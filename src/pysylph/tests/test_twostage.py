import io
import os
import pathlib
import random
import shutil
import tempfile
import unittest
import importlib.resources

from pysylph import (
    Database,
    GenomeSketch,
    Profiler,
    SampleSketch,
    Sketcher,
    TwoStageDatabase,
)


def data_path(name):
    # the fixtures are plain files next to the tests
    if hasattr(importlib.resources, "files"):
        return os.fspath(importlib.resources.files(__package__).joinpath(name))
    return os.fspath(pathlib.Path(__file__).parent / name)


def read_bytes(name):
    with open(data_path(name), "rb") as f:
        return f.read()


def random_genome(size, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(size))


class _TwoStageTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db = Database.load(data_path("ecoli.syldb"))
        # `ecoli.syl2db` was generated with the `sylph` v1.0.0 binary:
        #   sylph convert-db-two-screen ecoli.syldb -o ecoli
        cls.syl2db = read_bytes("ecoli.syl2db")
        cls.tsdb = TwoStageDatabase.load(data_path("ecoli.syl2db"))

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def path(self, name):
        return os.path.join(self.tmp, name)


class TestTwoStageDatabase(_TwoStageTestCase):
    def test_properties(self):
        self.assertEqual(len(self.tsdb), 3)
        self.assertEqual(self.tsdb.c, 200)
        self.assertEqual(self.tsdb.k, 31)
        self.assertEqual(self.tsdb.screen_c, 3000)

    def test_getitem(self):
        for i in range(3):
            sketch = self.tsdb[i]
            expected = self.db[i]
            self.assertIsInstance(sketch, GenomeSketch)
            self.assertEqual(sketch.name, expected.name)
            self.assertEqual(sketch.description, expected.description)
            self.assertEqual(sketch.genome_size, expected.genome_size)
            self.assertEqual(sketch.min_spacing, expected.min_spacing)
            self.assertEqual(sketch.c, 200)
            self.assertEqual(sketch.k, 31)
            self.assertEqual(sketch.kmers, sorted(expected.kmers))

    def test_getitem_negative(self):
        self.assertEqual(self.tsdb[-1].name, self.tsdb[2].name)
        self.assertEqual(self.tsdb[-3].name, self.tsdb[0].name)

    def test_getitem_index_error(self):
        with self.assertRaises(IndexError):
            self.tsdb[3]
        with self.assertRaises(IndexError):
            self.tsdb[-4]

    def test_iter(self):
        self.assertEqual(
            [g.name for g in self.tsdb], [g.name for g in self.db]
        )

    def test_load_path(self):
        db = TwoStageDatabase.load(pathlib.Path(data_path("ecoli.syl2db")))
        self.assertEqual(len(db), 3)

    def test_load_file_object(self):
        with open(data_path("ecoli.syl2db"), "rb") as f:
            db = TwoStageDatabase.load(f)
        self.assertEqual(len(db), 3)
        self.assertEqual(db[1].kmers, self.tsdb[1].kmers)
        db = TwoStageDatabase.load(io.BytesIO(self.syl2db))
        self.assertEqual(len(db), 3)

    def test_load_errors(self):
        with self.assertRaises(OSError):
            TwoStageDatabase.load(self.path("missing.syl2db"))
        garbage = b"NOPE" + bytes(range(256)) * 4
        for data in (b"", b"SY2D", garbage, self.syl2db[:100], self.syl2db[:-10]):
            with self.subTest(size=len(data)):
                with self.assertRaises(ValueError):
                    TwoStageDatabase.load(io.BytesIO(data))
                path = self.path("bad.syl2db")
                with open(path, "wb") as f:
                    f.write(data)
                with self.assertRaises(ValueError):
                    TwoStageDatabase.load(path)

    def test_load_syldb_is_error(self):
        with self.assertRaises(ValueError):
            TwoStageDatabase.load(data_path("ecoli.syldb"))

    def test_from_database_identical_to_cli(self):
        db = TwoStageDatabase.from_database(self.db)
        path = self.path("py.syl2db")
        db.dump(path)
        with open(path, "rb") as f:
            self.assertEqual(f.read(), self.syl2db)
        self.assertEqual(len(db), 3)
        self.assertEqual(db.c, self.tsdb.c)
        self.assertEqual(db.k, self.tsdb.k)
        self.assertEqual(db.screen_c, self.tsdb.screen_c)

    def test_from_database_screen_c(self):
        db = TwoStageDatabase.from_database(self.db, screen_c=1000)
        self.assertEqual(db.screen_c, 1000)
        self.assertEqual(len(db), 3)

    def test_dump_load_roundtrip(self):
        # dump of a database loaded from a path (copy of the file)
        path = self.path("copy.syl2db")
        self.tsdb.dump(path)
        with open(path, "rb") as f:
            self.assertEqual(f.read(), self.syl2db)
        # dump of a database loaded in memory, to a `PathLike`
        with open(path, "rb") as f:
            db = TwoStageDatabase.load(f)
        path2 = pathlib.Path(self.path("copy2.syl2db"))
        db.dump(path2)
        self.assertEqual(path2.read_bytes(), self.syl2db)
        # and load them back from a path and a file object
        reloaded = TwoStageDatabase.load(path2)
        self.assertEqual(reloaded[2].kmers, self.tsdb[2].kmers)
        with open(path2, "rb") as f:
            reloaded = TwoStageDatabase.load(f)
        self.assertEqual(reloaded[2].kmers, self.tsdb[2].kmers)

    def test_dump_survives_source_removal(self):
        path = self.path("source.syl2db")
        shutil.copy(data_path("ecoli.syl2db"), path)
        db = TwoStageDatabase.load(path)
        os.remove(path)
        db.dump(self.path("dest.syl2db"))
        with open(self.path("dest.syl2db"), "rb") as f:
            self.assertEqual(f.read(), self.syl2db)

    def test_dump_over_source(self):
        path = self.path("source.syl2db")
        shutil.copy(data_path("ecoli.syl2db"), path)
        db = TwoStageDatabase.load(path)
        db.dump(path)
        with open(path, "rb") as f:
            self.assertEqual(f.read(), self.syl2db)
        self.assertEqual(len(db[0].kmers), len(self.tsdb[0].kmers))


class TestTwoStageFromDatabase(_TwoStageTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.sketcher = Sketcher(c=20)
        cls.genome1 = cls.sketcher.sketch_genome("g1", [random_genome(50000, 1)])
        cls.genome2 = cls.sketcher.sketch_genome("g2", [random_genome(50000, 2)])

    def test_small_database(self):
        db = TwoStageDatabase.from_database(
            Database([self.genome1, self.genome2]), screen_c=100
        )
        self.assertEqual(len(db), 2)
        self.assertEqual(db.c, 20)
        self.assertEqual(db.screen_c, 100)
        self.assertEqual(db[1].name, "g2")
        self.assertEqual(db[1].kmers, sorted(self.genome2.kmers))

    def test_min_contain_drops_genomes(self):
        tiny = self.sketcher.sketch_genome("tiny", ["ACGT" * 5])
        database = Database([self.genome1, tiny, self.genome2])
        self.assertEqual(len(TwoStageDatabase.from_database(database)), 2)
        db = TwoStageDatabase.from_database(database, min_contain=0)
        self.assertEqual(len(db), 3)

    def test_empty(self):
        with self.assertRaises(ValueError):
            TwoStageDatabase.from_database(Database())

    def test_all_dropped(self):
        with self.assertRaises(ValueError):
            TwoStageDatabase.from_database(
                Database([self.genome1]), min_contain=1000000
            )

    def test_inconsistent_c(self):
        other = Sketcher(c=40).sketch_genome("g3", [random_genome(50000, 3)])
        with self.assertRaises(ValueError):
            TwoStageDatabase.from_database(Database([self.genome1, other]))

    def test_inconsistent_k(self):
        other = Sketcher(c=20, k=21).sketch_genome("g3", [random_genome(50000, 3)])
        with self.assertRaises(ValueError):
            TwoStageDatabase.from_database(Database([self.genome1, other]))

    def test_no_profiling(self):
        other = self.sketcher.sketch_genome(
            "g3", [random_genome(50000, 3)], profiling=False
        )
        with self.assertRaises(ValueError):
            TwoStageDatabase.from_database(Database([self.genome1, other]))

    def test_screen_c_too_small(self):
        database = Database([self.genome1])
        with self.assertRaises(ValueError):
            TwoStageDatabase.from_database(database, screen_c=10)
        TwoStageDatabase.from_database(database, screen_c=20)

    def test_min_sparse_kmers(self):
        with self.assertRaises(ValueError):
            TwoStageDatabase.from_database(Database([self.genome1]), min_sparse_kmers=0)

    def test_type_error(self):
        with self.assertRaises(TypeError):
            TwoStageDatabase.from_database([self.genome1])

    def test_profile_built_database(self):
        # sample made of reads from the second genome
        genome = random_genome(50000, 2)
        reads = [genome[i : i + 150] for i in range(0, 49000, 50)]
        sample = self.sketcher.sketch_single("reads", reads)
        db = TwoStageDatabase.from_database(
            Database([self.genome1, self.genome2]), screen_c=100
        )
        results = Profiler(db).query(sample)
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].genome_sketch.name, "g2")
        self.assertEqual(results[0].genome_sketch.kmers, sorted(self.genome2.kmers))
        results = Profiler(db).profile(sample)
        self.assertEqual(len(results), 1)
        self.assertAlmostEqual(results[0].sequence_abundance, 100.0)


class TestTwoStageProfiler(_TwoStageTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.o157 = SampleSketch.load(data_path("o157.sylsp"))
        cls.k12 = SampleSketch.load(data_path("k12.paired.sylsp"))

    # Expected values are the ones of the `sylph` v1.0.0 binary:
    #   sylph query ecoli.syl2db o157.sylsp
    #   sylph profile ecoli.syl2db o157.sylsp
    # which are the same as with `ecoli.syldb`.

    def test_query(self):
        results = Profiler(self.tsdb).query(self.o157)
        self.assertEqual(len(results), 3)
        for result, name, ani, naive, cov in zip(
            results,
            ["e.coli-o157", "e.coli-EC590", "e.coli-K12"],
            [99.79, 98.32, 98.22],
            [96.02, 94.29, 94.26],
            [0.360, 0.319, 0.327],
        ):
            self.assertEqual(result.genome_sketch.name, f"test_files/{name}.fasta.gz")
            self.assertAlmostEqual(result.ani, ani, places=2)
            self.assertAlmostEqual(result.ani_naive, naive, places=2)
            self.assertAlmostEqual(result.coverage, cov, places=3)

    def test_query_paired(self):
        results = Profiler(self.tsdb).query(self.k12)
        self.assertEqual(len(results), 3)
        for result, name, ani, naive, cov in zip(
            results,
            ["K12", "EC590", "o157"],
            [98.89, 98.11, 97.09],
            [88.48, 88.15, 86.73],
            [0.032, 0.037, 0.031],
        ):
            self.assertEqual(result.genome_sketch.name, f"test_files/e.coli-{name}.fasta.gz")
            self.assertAlmostEqual(result.ani, ani, places=2)
            self.assertAlmostEqual(result.ani_naive, naive, places=2)
            self.assertAlmostEqual(result.coverage, cov, places=3)

    def test_profile(self):
        results = Profiler(self.tsdb).profile(self.o157)
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].genome_sketch.name, "test_files/e.coli-o157.fasta.gz")
        self.assertAlmostEqual(results[0].ani, 99.79, places=2)
        self.assertAlmostEqual(results[0].ani_naive, 96.02, places=2)
        self.assertAlmostEqual(results[0].coverage, 0.360, places=3)
        self.assertAlmostEqual(results[0].sequence_abundance, 100.0, places=2)
        self.assertAlmostEqual(results[0].taxonomic_abundance, 100.0, places=2)
        self.assertEqual(results[0].kmers_reassigned, 0)
        self.assertAlmostEqual(results[0].ani_ci[0], 99.64, places=2)
        self.assertAlmostEqual(results[0].ani_ci[1], 99.95, places=2)
        self.assertAlmostEqual(results[0].effective_lambda, 0.360, places=3)
        self.assertAlmostEqual(results[0].lambda_ci[0], 0.34, places=2)
        self.assertAlmostEqual(results[0].lambda_ci[1], 0.38, places=2)
        self.assertAlmostEqual(results[0].mean_coverage, 1.187, places=3)
        self.assertEqual(results[0].containment, (6208, 21899))

    def test_profile_estimate_unknown(self):
        # sylph profile ecoli.syl2db k12.paired.sylsp -u
        profiler = Profiler(self.tsdb, estimate_unknown=True)
        results = profiler.profile(self.k12)
        self.assertEqual(len(results), 1)
        self.assertAlmostEqual(results[0].ani, 98.89, places=2)
        self.assertAlmostEqual(results[0].coverage, 0.066, places=3)
        self.assertAlmostEqual(results[0].sequence_abundance, 87.7690, places=3)
        self.assertAlmostEqual(results[0].taxonomic_abundance, 100.0, places=3)

    def test_same_as_database(self):
        for sample in (self.o157, self.k12):
            for kwargs in ({}, {"estimate_unknown": True}, {"minimum_ani": 98.0}):
                with self.subTest(sample=sample.k, **kwargs):
                    plain = Profiler(self.db, **kwargs)
                    twostage = Profiler(self.tsdb, **kwargs)
                    for method in ("query", "profile"):
                        expected = getattr(plain, method)(sample)
                        actual = getattr(twostage, method)(sample)
                        self.assertEqual(len(actual), len(expected))
                        for a, e in zip(actual, expected):
                            self.assertEqual(a.genome_sketch.name, e.genome_sketch.name)
                            self.assertEqual(a.ani, e.ani)
                            self.assertEqual(a.ani_naive, e.ani_naive)
                            self.assertEqual(a.coverage, e.coverage)
                            if method == "profile":
                                self.assertEqual(a.sequence_abundance, e.sequence_abundance)
                                self.assertEqual(a.taxonomic_abundance, e.taxonomic_abundance)
                                self.assertEqual(a.kmers_reassigned, e.kmers_reassigned)

    def test_in_memory_database(self):
        db = TwoStageDatabase.from_database(self.db)
        results = Profiler(db).query(self.o157)
        self.assertEqual(len(results), 3)
        self.assertAlmostEqual(results[0].ani, 99.79, places=2)

    def test_database_property(self):
        self.assertIs(Profiler(self.tsdb).database, self.tsdb)
        self.assertIs(Profiler(self.db).database, self.db)

    def test_screen_ani(self):
        # the o157 genome has an adjusted ANI of 99.79, the others ~98.3
        profiler = Profiler(self.tsdb, screen_ani=99.0)
        self.assertEqual(len(profiler.query(self.o157)), 1)
        profiler = Profiler(self.tsdb, screen_ani=100.0)
        self.assertEqual(profiler.query(self.o157), [])
        self.assertEqual(profiler.profile(self.o157), [])
        # ignored for plain databases
        profiler = Profiler(self.db, screen_ani=100.0)
        self.assertEqual(len(profiler.query(self.o157)), 3)

    def test_screen_ani_invalid(self):
        with self.assertRaises(ValueError):
            Profiler(self.tsdb, screen_ani=-1.0)
        with self.assertRaises(ValueError):
            Profiler(self.tsdb, screen_ani=101.0)

    def test_min_contain(self):
        profiler = Profiler(self.tsdb, min_contain=3120)
        self.assertEqual(len(profiler.query(self.o157)), 2)
        profiler = Profiler(self.tsdb, min_contain=10000)
        self.assertEqual(profiler.query(self.o157), [])

    def test_k_mismatch(self):
        sample = Sketcher(k=21).sketch_single("reads", ["ACGT" * 100])
        profiler = Profiler(self.tsdb)
        with self.assertRaises(ValueError):
            profiler.query(sample)
        with self.assertRaises(ValueError):
            profiler.profile(sample)

    def test_c_mismatch(self):
        sample = Sketcher(c=1000).sketch_single("reads", ["ACGT" * 100])
        profiler = Profiler(self.tsdb)
        with self.assertRaises(ValueError):
            profiler.query(sample)
        with self.assertRaises(ValueError):
            profiler.profile(sample)

    def test_invalid_database(self):
        with self.assertRaises(TypeError):
            Profiler([self.db[0]])
        with self.assertRaises(TypeError):
            Profiler(None)
