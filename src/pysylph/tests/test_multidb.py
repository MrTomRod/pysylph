import gzip
import os
import pathlib
import unittest
import importlib.resources
import warnings

from pysylph import (
    Database,
    GenomeSketch,
    Profiler,
    SampleSketch,
    Sketcher,
    TwoStageDatabase,
)


def data_path(name):
    if hasattr(importlib.resources, "files"):
        return os.fspath(importlib.resources.files(__package__).joinpath(name))
    return os.fspath(pathlib.Path(__file__).parent / name)


class TestMultiDatabaseProfiler(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ecoli = Database.load(data_path("ecoli.syldb"))
        cls.clodf13 = Database.load(data_path("clodf13.syldb"))
        cls.ecoli2 = TwoStageDatabase.load(data_path("ecoli.syl2db"))
        cls.o157 = SampleSketch.load(data_path("o157.sylsp"))
        cls.k12 = SampleSketch.load(data_path("k12.paired.sylsp"))
        # K-12 reads (~0.03x) mixed with CloDF13 plasmid reads (~3x)
        cls.mix = SampleSketch.load(data_path("mix.sylsp"))
        with gzip.open(data_path("clodf13.fna.gz"), "rt") as f:
            cls.clodf13_seq = "".join(l.strip() for l in f if not l.startswith(">"))
        cls.merged = Database(list(cls.ecoli) + list(cls.clodf13))

    def assertSameResults(self, actual, expected):
        self.assertEqual(len(actual), len(expected))
        for a, e in zip(actual, expected):
            self.assertEqual(a.genome_sketch.name, e.genome_sketch.name)
            self.assertEqual(a.ani, e.ani)
            self.assertEqual(a.ani_naive, e.ani_naive)
            self.assertEqual(a.coverage, e.coverage)
            if hasattr(e, "sequence_abundance"):
                self.assertEqual(a.sequence_abundance, e.sequence_abundance)
                self.assertEqual(a.taxonomic_abundance, e.taxonomic_abundance)
                self.assertEqual(a.kmers_reassigned, e.kmers_reassigned)

    # Expected values are the ones of the `sylph` binary (expose-library-api):
    #   sylph query   ecoli.syl2db clodf13.syldb mix.sylsp
    #   sylph profile ecoli.syl2db clodf13.syldb mix.sylsp [-u]
    # which are the same as with `ecoli.syldb clodf13.syldb`, in any order.
    # mix.sylsp was made with `sylph sketch -r` of the reads of k12.R1/R2.fq.gz
    # and reads randomly sampled from the clodf13 genome.

    def test_properties(self):
        profiler = Profiler([self.ecoli, self.clodf13])
        self.assertEqual(profiler.databases, (self.ecoli, self.clodf13))
        self.assertIsInstance(profiler.databases, tuple)
        profiler = Profiler((self.ecoli2, self.clodf13))
        self.assertIs(profiler.databases[0], self.ecoli2)
        self.assertIs(profiler.databases[1], self.clodf13)

    def test_database_property(self):
        dbs = [self.ecoli2, self.clodf13]
        self.assertIs(Profiler(dbs).database, dbs)
        dbs = (self.ecoli2, self.clodf13)
        self.assertIs(Profiler(dbs).database, dbs)
        self.assertIs(Profiler(self.ecoli).database, self.ecoli)
        self.assertEqual(Profiler(self.ecoli).databases, (self.ecoli,))
        self.assertIs(Profiler(self.ecoli2).database, self.ecoli2)
        self.assertEqual(Profiler(self.ecoli2).databases, (self.ecoli2,))

    def test_same_as_merged_database(self):
        for i, sample in enumerate((self.o157, self.k12, self.mix)):
            for kwargs in ({}, {"estimate_unknown": True}, {"minimum_ani": 98.0}):
                with self.subTest(sample=i, **kwargs):
                    expected = Profiler(self.merged, **kwargs)
                    actual = Profiler([self.ecoli, self.clodf13], **kwargs)
                    reverse = Profiler([self.clodf13, self.ecoli], **kwargs)
                    for method in ("query", "profile"):
                        e = getattr(expected, method)(sample)
                        self.assertSameResults(getattr(actual, method)(sample), e)
                        self.assertEqual(
                            sorted(r.genome_sketch.name for r in e),
                            sorted(r.genome_sketch.name for r in getattr(reverse, method)(sample)),
                        )

    def test_same_as_merged_database_twostage(self):
        merged2 = TwoStageDatabase.from_database(self.merged)
        for i, sample in enumerate((self.o157, self.k12, self.mix)):
            with self.subTest(sample=i):
                a = Profiler([self.ecoli2, self.clodf13])
                b = Profiler(merged2)
                for method in ("query", "profile"):
                    self.assertSameResults(
                        getattr(a, method)(sample), getattr(b, method)(sample)
                    )

    def test_two_twostage_databases(self):
        clodf13_2 = TwoStageDatabase.from_database(self.clodf13)
        profiler = Profiler([self.ecoli2, clodf13_2])
        for method in ("query", "profile"):
            self.assertSameResults(
                getattr(profiler, method)(self.mix),
                getattr(Profiler(self.merged), method)(self.mix),
            )

    def test_mixed_query_mix(self):
        for dbs in (
            [self.ecoli2, self.clodf13],
            [self.clodf13, self.ecoli2],
            [self.ecoli, self.clodf13],
        ):
            results = Profiler(dbs).query(self.mix)
            self.assertEqual(len(results), 4)
            for result, name, ani, naive, cov in zip(
                results,
                [
                    "src/pysylph/tests/clodf13.fna.gz",
                    "test_files/e.coli-K12.fasta.gz",
                    "test_files/e.coli-EC590.fasta.gz",
                    "test_files/e.coli-o157.fasta.gz",
                ],
                [100.00, 98.89, 98.11, 97.09],
                [99.69, 88.48, 88.15, 86.73],
                [2.200, 0.032, 0.037, 0.031],
            ):
                self.assertEqual(result.genome_sketch.name, name)
                self.assertAlmostEqual(result.ani, ani, places=2)
                self.assertAlmostEqual(result.ani_naive, naive, places=2)
                self.assertAlmostEqual(result.coverage, cov, places=3)

    def test_mixed_profile_mix(self):
        for dbs in (
            [self.ecoli2, self.clodf13],
            [self.clodf13, self.ecoli2],
            [self.ecoli, self.clodf13],
        ):
            results = Profiler(dbs).profile(self.mix)
            self.assertEqual(len(results), 2)
            for result, name, ani, naive, cov, tax, seq in zip(
                results,
                ["src/pysylph/tests/clodf13.fna.gz", "test_files/e.coli-K12.fasta.gz"],
                [100.00, 98.89],
                [99.69, 88.48],
                [2.200, 0.032],
                [98.5483, 1.4517],
                [12.7002, 87.2998],
            ):
                self.assertEqual(result.genome_sketch.name, name)
                self.assertAlmostEqual(result.ani, ani, places=2)
                self.assertAlmostEqual(result.ani_naive, naive, places=2)
                self.assertAlmostEqual(result.coverage, cov, places=3)
                self.assertAlmostEqual(result.taxonomic_abundance, tax, places=3)
                self.assertAlmostEqual(result.sequence_abundance, seq, places=3)
                self.assertEqual(result.kmers_reassigned, 0)

    def test_mixed_profile_estimate_unknown(self):
        profiler = Profiler([self.ecoli2, self.clodf13], estimate_unknown=True)
        results = profiler.profile(self.mix)
        self.assertEqual(len(results), 2)
        self.assertAlmostEqual(results[0].taxonomic_abundance, 98.5483, places=3)
        self.assertAlmostEqual(results[0].sequence_abundance, 11.7054, places=3)
        self.assertAlmostEqual(results[1].taxonomic_abundance, 1.4517, places=3)
        self.assertAlmostEqual(results[1].sequence_abundance, 80.4616, places=3)

    def test_mixed_o157_k12(self):
        # same as `sylph query|profile ecoli.syl2db clodf13.syldb sample`
        for dbs in ([self.ecoli2], [self.ecoli2, self.clodf13], (self.clodf13, self.ecoli2)):
            profiler = Profiler(dbs)
            results = profiler.query(self.o157)
            self.assertEqual(len(results), 3)
            for result, name, ani, naive, cov in zip(
                results,
                ["o157", "EC590", "K12"],
                [99.79, 98.32, 98.22],
                [96.02, 94.29, 94.26],
                [0.360, 0.319, 0.327],
            ):
                self.assertEqual(result.genome_sketch.name, f"test_files/e.coli-{name}.fasta.gz")
                self.assertAlmostEqual(result.ani, ani, places=2)
                self.assertAlmostEqual(result.ani_naive, naive, places=2)
                self.assertAlmostEqual(result.coverage, cov, places=3)
            results = profiler.profile(self.o157)
            self.assertEqual(len(results), 1)
            self.assertEqual(results[0].genome_sketch.name, "test_files/e.coli-o157.fasta.gz")
            self.assertAlmostEqual(results[0].ani, 99.79, places=2)
            self.assertAlmostEqual(results[0].sequence_abundance, 100.0, places=2)

            results = profiler.query(self.k12)
            self.assertEqual(
                [r.genome_sketch.name for r in results],
                [f"test_files/e.coli-{n}.fasta.gz" for n in ("K12", "EC590", "o157")],
            )
            for result, ani, naive in zip(
                results, [98.89, 98.11, 97.09], [88.48, 88.15, 86.73]
            ):
                self.assertAlmostEqual(result.ani, ani, places=2)
                self.assertAlmostEqual(result.ani_naive, naive, places=2)

    def test_single_item_sequence(self):
        for db in (self.ecoli, self.ecoli2):
            for method in ("query", "profile"):
                self.assertSameResults(
                    getattr(Profiler([db]), method)(self.o157),
                    getattr(Profiler(db), method)(self.o157),
                )

    def test_no_warning(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            Profiler([self.ecoli2, self.clodf13])
            Profiler(self.ecoli)

    def test_empty(self):
        with self.assertRaises(ValueError):
            Profiler([])
        with self.assertRaises(ValueError):
            Profiler(())

    def test_invalid_database(self):
        with self.assertRaises(TypeError):
            Profiler([self.ecoli, None])
        with self.assertRaises(TypeError):
            Profiler([self.ecoli, self.ecoli[0]])
        with self.assertRaises(TypeError):
            Profiler([[self.ecoli]])
        with self.assertRaises(TypeError):
            Profiler("ecoli.syldb")
        with self.assertRaises(TypeError):
            Profiler(iter([self.ecoli]))

    def test_k_mismatch(self):
        sketch = Sketcher(k=21).sketch_genome("clodf13_k21", [self.clodf13_seq])
        self.assertEqual(sketch.k, 21)
        db21 = Database([sketch])
        for dbs in (
            [self.ecoli, db21],
            [self.ecoli2, db21],
            [db21, self.ecoli2],
            (self.ecoli2, self.clodf13, db21),
        ):
            with self.assertRaises(ValueError):
                Profiler(dbs)
        # same k is fine
        Profiler([db21])

    def test_c_not_required_to_match(self):
        sequence = self.clodf13_seq
        db100 = Database([Sketcher(c=100).sketch_genome("clodf13_c100", [sequence])])
        db300 = Database([Sketcher(c=300).sketch_genome("clodf13_c300", [sequence])])
        # the sample c must not exceed the minimum c across all databases
        for dbs, ok, too_high in (
            ([self.ecoli2, db100], 100, 150),
            ([db100, self.ecoli2], 100, 150),
            ([self.ecoli, db100], 100, 150),
            ([self.ecoli2, db300], 200, 250),
            ([db300, self.ecoli], 200, 250),
        ):
            profiler = Profiler(dbs)
            sample = Sketcher(c=ok).sketch_single("reads", [sequence])
            profiler.query(sample)
            profiler.profile(sample)
            sample = Sketcher(c=too_high).sketch_single("reads", [sequence])
            with self.assertRaises(ValueError):
                profiler.query(sample)
            with self.assertRaises(ValueError):
                profiler.profile(sample)

    def test_sample_checks_all_databases(self):
        sample = Sketcher(c=1000).sketch_single("reads", ["ACGT" * 100])
        for dbs in ([self.ecoli, self.clodf13], [self.ecoli2, self.clodf13], [self.clodf13, self.ecoli2]):
            profiler = Profiler(dbs)
            with self.assertRaises(ValueError):
                profiler.query(sample)
            with self.assertRaises(ValueError):
                profiler.profile(sample)
        sample = Sketcher(k=21).sketch_single("reads", ["ACGT" * 100])
        with self.assertRaises(ValueError):
            Profiler([self.ecoli2, self.clodf13]).query(sample)

    def test_duplicate_genomes(self):
        for dbs in (
            [self.ecoli, self.ecoli2],
            [self.ecoli2, self.ecoli],
            [self.ecoli2, self.clodf13, self.ecoli2],
            [self.ecoli, self.ecoli],
        ):
            with self.assertWarns(UserWarning) as cm:
                Profiler(dbs)
            self.assertIn("more than once", str(cm.warning))
            self.assertIn("e.coli-o157", str(cm.warning))
        # warns once, at construction
        with self.assertWarns(UserWarning):
            profiler = Profiler([self.ecoli, self.ecoli2])
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            profiler.query(self.o157)


if __name__ == "__main__":
    unittest.main()
