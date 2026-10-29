import os
import pathlib
import unittest
import gzip
import random
import importlib.resources
from contextlib import nullcontext

from pysylph import Database, SampleSketch, Sketcher, Profiler


def data_path(name):
    if hasattr(importlib.resources, "files"):
        return os.fspath(importlib.resources.files(__package__).joinpath(name))
    return os.fspath(pathlib.Path(__file__).parent / name)


class TestProfiler(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # load reference database
        if hasattr(importlib.resources, "files"):
            handler = nullcontext(
                importlib.resources.files(__package__).joinpath("ecoli.syldb")
            )
        else:
            handler = importlib.resources.path(__package__, "ecoli.syldb")
        with handler as path:
            cls.db = Database.load(path)
        # load seq
        if hasattr(importlib.resources, "files"):
            handler = nullcontext(
                importlib.resources.files(__package__).joinpath("o157.sylsp")
            )
        else:
            handler = importlib.resources.path(__package__, "o157.sylsp")
        with handler as path:
            cls.sample = SampleSketch.load(path)

    def test_query(self):
        # sketch sequence
        profiler = Profiler(self.db)
        results = profiler.query(self.sample)

        self.assertEqual(len(results), 3)
        self.assertAlmostEqual(results[0].ani, 99.79, places=2)
        self.assertAlmostEqual(results[0].ani_naive, 96.02, places=2)

        self.assertAlmostEqual(results[1].ani, 98.32, places=2)
        self.assertAlmostEqual(results[1].ani_naive, 94.29, places=2)

        self.assertAlmostEqual(results[2].ani, 98.22, places=2)
        self.assertAlmostEqual(results[2].ani_naive, 94.26, places=2)

    def test_profile(self):
        # sketch sequence
        profiler = Profiler(self.db)
        results = profiler.profile(self.sample)

        self.assertEqual(len(results), 1)
        self.assertAlmostEqual(results[0].ani, 99.79, places=2)
        self.assertAlmostEqual(results[0].ani_naive, 96.02, places=2)
        self.assertAlmostEqual(results[0].sequence_abundance, 100.0, places=2)
        self.assertAlmostEqual(results[0].taxonomic_abundance, 100.0, places=2)
        self.assertAlmostEqual(results[0].kmers_reassigned, 0)

    def test_fields_query(self):
        # expected values from `sylph query ecoli.syldb o157.sylsp`
        results = Profiler(self.db).query(self.sample)
        expected = [
            # contig, ani_ci, lambda, lambda_ci, mean_cov, containment
            ("NZ_CP017438.1 Escherichia coli O157:H7 strain 2159 chromosome, complete genome",
             (99.64, 99.95), 0.360, (0.34, 0.38), 1.187, (6208, 21899)),
            ("NZ_CP016182.2 Escherichia coli strain EC590 chromosome, complete genome",
             (98.05, 98.58), 0.319, (0.29, 0.35), 1.172, (3122, 19330)),
            ("NC_007779.1 Escherichia coli str. K-12 substr. W3110, complete sequence",
             (97.91, 98.47), 0.327, (0.30, 0.36), 1.171, (3114, 19485)),
        ]
        self.assertEqual(len(results), len(expected))
        for result, (contig, ani_ci, lam, lam_ci, mean_cov, containment) in zip(results, expected):
            self.assertEqual(result.genome_sketch.description, contig)
            self.assertEqual(len(result.ani_ci), 2)
            self.assertAlmostEqual(result.ani_ci[0], ani_ci[0], places=2)
            self.assertAlmostEqual(result.ani_ci[1], ani_ci[1], places=2)
            self.assertEqual(result.lambda_status, "OK")
            self.assertAlmostEqual(result.effective_lambda, lam, places=3)
            self.assertAlmostEqual(result.lambda_ci[0], lam_ci[0], places=2)
            self.assertAlmostEqual(result.lambda_ci[1], lam_ci[1], places=2)
            self.assertAlmostEqual(result.median_coverage, 1.0)
            self.assertAlmostEqual(result.mean_coverage, mean_cov, places=3)
            self.assertEqual(result.containment, containment)
            # without `-u`, the coverage is the effective lambda
            self.assertAlmostEqual(result.coverage, lam, places=3)

    def test_fields_profile(self):
        # expected values from `sylph profile ecoli.syldb o157.sylsp`
        (result,) = Profiler(self.db).profile(self.sample)
        self.assertEqual(
            result.genome_sketch.description,
            "NZ_CP017438.1 Escherichia coli O157:H7 strain 2159 chromosome, complete genome",
        )
        self.assertAlmostEqual(result.ani_ci[0], 99.64, places=2)
        self.assertAlmostEqual(result.ani_ci[1], 99.95, places=2)
        self.assertEqual(result.lambda_status, "OK")
        self.assertAlmostEqual(result.effective_lambda, 0.360, places=3)
        self.assertAlmostEqual(result.lambda_ci[0], 0.34, places=2)
        self.assertAlmostEqual(result.lambda_ci[1], 0.38, places=2)
        self.assertAlmostEqual(result.median_coverage, 1.0)
        self.assertAlmostEqual(result.mean_coverage, 1.187, places=3)
        self.assertEqual(result.containment, (6208, 21899))

    def test_fields_estimate_unknown(self):
        # expected values from `sylph query|profile -u ecoli.syldb o157.sylsp`:
        # `Eff_cov` becomes `True_cov`, but lambda is unchanged
        profiler = Profiler(self.db, estimate_unknown=True)
        for results in (profiler.query(self.sample), profiler.profile(self.sample)):
            self.assertAlmostEqual(results[0].coverage, 3.108, places=3)
            self.assertAlmostEqual(results[0].effective_lambda, 0.360, places=3)
            self.assertAlmostEqual(results[0].ani_ci[0], 99.64, places=2)
            self.assertEqual(results[0].containment, (6208, 21899))
        results = profiler.query(self.sample)
        self.assertAlmostEqual(results[1].coverage, 2.758, places=3)
        self.assertAlmostEqual(results[1].effective_lambda, 0.319, places=3)
        self.assertAlmostEqual(results[2].coverage, 2.821, places=3)
        self.assertAlmostEqual(results[2].effective_lambda, 0.327, places=3)

    def test_fields_ci_unclamped(self):
        # `ani` is capped at 100 but the CI is printed as-is by `sylph`:
        # `sylph query ecoli.syldb k12.paired.sylsp` gives 96.80-100.69
        sample = SampleSketch.load(
            importlib.resources.files(__package__).joinpath("k12.paired.sylsp")
        )
        results = Profiler(self.db).query(sample)
        self.assertEqual(len(results), 3)
        self.assertAlmostEqual(results[0].ani, 98.89, places=2)
        self.assertAlmostEqual(results[0].ani_ci[0], 96.80, places=2)
        self.assertAlmostEqual(results[0].ani_ci[1], 100.69, places=2)
        self.assertAlmostEqual(results[0].lambda_ci[0], 0.02, places=2)
        self.assertAlmostEqual(results[0].lambda_ci[1], 0.06, places=2)
        self.assertAlmostEqual(results[0].mean_coverage, 1.016, places=3)
        self.assertEqual(results[0].containment, (439, 19485))
        self.assertEqual(results[1].containment, (387, 19330))
        self.assertEqual(results[2].containment, (265, 21899))

    def test_fields_low_lambda(self):
        # `sylph query -m 80 -M 1 --min-contain 1 ecoli.syldb k12.R1.sylsp`
        # reports `LOW` lambdas and no confidence intervals (`NA-NA`)
        sample = SampleSketch.load(
            importlib.resources.files(__package__).joinpath("k12.R1.sylsp")
        )
        profiler = Profiler(self.db, minimum_ani=80.0, min_number_kmers=1, min_contain=1)
        results = profiler.query(sample)
        self.assertEqual(len(results), 3)
        for result, ani, cov, mean_cov, containment in zip(
            results,
            [86.53, 86.29, 84.78],
            [1.009, 1.010, 1.008],
            [1.009, 1.010, 1.008],
            [(220, 19485), (200, 19330), (131, 21899)],
        ):
            self.assertAlmostEqual(result.ani, ani, places=2)
            self.assertEqual(result.lambda_status, "LOW")
            self.assertIsNone(result.effective_lambda)
            self.assertIsNone(result.ani_ci)
            self.assertIsNone(result.lambda_ci)
            self.assertAlmostEqual(result.coverage, cov, places=3)
            self.assertAlmostEqual(result.median_coverage, 1.0)
            self.assertAlmostEqual(result.mean_coverage, mean_cov, places=3)
            self.assertEqual(result.containment, containment)

    def test_seq_id(self):
        # expected values from `sylph query|profile -u -I 99.5 ecoli.syldb o157.sylsp`
        profiler = Profiler(self.db, estimate_unknown=True, seq_id=99.5)
        results = profiler.query(self.sample)
        self.assertEqual(len(results), 3)
        self.assertAlmostEqual(results[0].ani, 99.79, places=2)
        self.assertAlmostEqual(results[0].coverage, 0.421, places=3)
        self.assertAlmostEqual(results[1].coverage, 0.374, places=3)
        self.assertAlmostEqual(results[2].coverage, 0.382, places=3)
        results = profiler.profile(self.sample)
        self.assertEqual(len(results), 1)
        self.assertAlmostEqual(results[0].coverage, 0.421, places=3)
        self.assertAlmostEqual(results[0].taxonomic_abundance, 100.0, places=4)
        self.assertAlmostEqual(results[0].sequence_abundance, 42.2897, places=4)

    def test_min_contain(self):
        # o157 hit has 6208 contained k-mers, others have 3122 and 3114
        profiler = Profiler(self.db, min_contain=3120)
        results = profiler.query(self.sample)
        self.assertEqual(len(results), 2)
        profiler = Profiler(self.db, min_contain=10000)
        self.assertEqual(profiler.query(self.sample), [])
        self.assertEqual(profiler.profile(self.sample), [])

    def test_min_number_kmers(self):
        profiler = Profiler(self.db, min_number_kmers=10000000)
        self.assertEqual(profiler.query(self.sample), [])

    def test_k_mismatch(self):
        sketcher = Sketcher(k=21)
        sample = sketcher.sketch_single("reads", ["ACGT" * 100])
        profiler = Profiler(self.db)
        with self.assertRaises(ValueError):
            profiler.query(sample)
        with self.assertRaises(ValueError):
            profiler.profile(sample)

    def test_c_mismatch(self):
        sketcher = Sketcher(c=1000)
        sample = sketcher.sketch_single("reads", ["ACGT" * 100])
        profiler = Profiler(self.db)
        with self.assertRaises(ValueError):
            profiler.query(sample)
        with self.assertRaises(ValueError):
            profiler.profile(sample)

    def test_same_name(self):
        # genomes sharing a name must not be mixed up in the results
        rng = random.Random(0)
        genomes = ["".join(rng.choice("ACGT") for _ in range(50000)) for _ in range(2)]
        sketcher = Sketcher(c=20)
        sketches = [sketcher.sketch_genome("same", [g]) for g in genomes]
        reads = [genomes[1][i : i + 150] for i in range(0, 49000, 50)]
        sample = sketcher.sketch_single("reads", reads)
        profiler = Profiler(Database(sketches))
        for results in (profiler.query(sample), profiler.profile(sample)):
            self.assertEqual(len(results), 1)
            self.assertEqual(results[0].genome_sketch.kmers, sketches[1].kmers)


class TestEstimateReadCounts(unittest.TestCase):
    """``Profiler(estimate_read_counts=True)`` reports estimated read counts.

    This exposes the ``sylph --estimate-read-counts`` option: instead of a
    relative abundance, ``ProfileResult.sequence_abundance`` holds the
    estimated number of reads assigned to each genome, computed as

        True_cov * genome_size / mean_read_length * bases_explained

    and rounded, as in ``sylph``. Like ``sylph``, the option implies
    ``estimate_unknown`` (``-u``).
    """

    @classmethod
    def setUpClass(cls):
        cls.sketcher = Sketcher(c=20)
        cls.db = Database.load(data_path("ecoli.syldb"))
        cls.clodf13 = Database.load(data_path("clodf13.syldb"))
        cls.sample = SampleSketch.load(data_path("o157.sylsp"))
        # K-12 reads (~0.03x) mixed with CloDF13 plasmid reads (~3x)
        cls.mix = SampleSketch.load(data_path("mix.sylsp"))

    def test_profile_reference_values(self):
        # expected values from
        # `sylph profile --estimate-read-counts ecoli.syldb clodf13.syldb mix.sylsp`
        profiler = Profiler(self.db, self.clodf13, estimate_read_counts=True)
        results = profiler.profile(self.mix)
        self.assertEqual(len(results), 2)
        for result, count, abundance in zip(
            results, [566.0, 3889.0], [98.5483, 1.4517]
        ):
            # counts are rounded, but still floats
            self.assertIsInstance(result.sequence_abundance, float)
            self.assertTrue(result.sequence_abundance.is_integer())
            self.assertAlmostEqual(result.sequence_abundance, count)
            # the taxonomic abundance is unchanged, and still a percentage
            self.assertAlmostEqual(result.taxonomic_abundance, abundance, places=4)
        # sorted by relative abundance, like `sylph`
        self.assertEqual(
            [r.genome_sketch.description.split(" ")[0] for r in results],
            ["X04466.1", "NC_007779.1"],
        )

    def test_single_hit(self):
        # `sylph profile --estimate-read-counts ecoli.syldb o157.sylsp`
        (result,) = Profiler(self.db, estimate_read_counts=True).profile(self.sample)
        self.assertAlmostEqual(result.sequence_abundance, 1149.0)
        self.assertAlmostEqual(result.taxonomic_abundance, 100.0)

    def test_implies_estimate_unknown(self):
        # `--estimate-read-counts` forces `-u`, as in `sylph`: the coverage
        # is the true coverage, not the effective coverage
        expected = Profiler(self.db, estimate_unknown=True)
        actual = Profiler(self.db, estimate_read_counts=True)
        # explicitly passing `estimate_unknown=False` is overridden
        explicit = Profiler(
            self.db, estimate_unknown=False, estimate_read_counts=True
        )
        for profiler in (actual, explicit):
            query = profiler.query(self.sample)
            self.assertEqual(
                [(r.coverage, r.ani) for r in query],
                [(r.coverage, r.ani) for r in expected.query(self.sample)],
            )
            # the coverage is the true coverage of `-u`, not the effective
            # coverage of the default
            self.assertAlmostEqual(query[0].coverage, 3.108, places=3)
            self.assertNotAlmostEqual(
                query[0].coverage,
                Profiler(self.db).query(self.sample)[0].coverage,
                places=3,
            )
            # the profile is filtered the same way as with `-u`
            self.assertEqual(
                [
                    (r.genome_sketch.name, r.taxonomic_abundance)
                    for r in profiler.profile(self.sample)
                ],
                [
                    (r.genome_sketch.name, r.taxonomic_abundance)
                    for r in expected.profile(self.sample)
                ],
            )

    def test_matches_formula(self):
        # check the counts against the formula used by `sylph`, with reads of
        # a known length (so that the mean read length is known exactly)
        rng = random.Random(3)
        genome = "".join(rng.choice("ACGT") for _ in range(100_000))
        database = Database([self.sketcher.sketch_genome("g", [genome])])
        read_length = 150
        reads = [
            genome[i : i + read_length]
            for i in range(0, len(genome) - read_length, read_length // 10)
        ]
        sample = self.sketcher.sketch_single("reads", reads)
        # `bases_explained` is the sum of the relative sequence abundances
        # scaled back by 100 (the abundances are normalized to 100%)
        (relative,) = Profiler(database, estimate_unknown=True).profile(sample)
        bases_explained = relative.sequence_abundance / 100.0
        (counts,) = Profiler(database, estimate_read_counts=True).profile(sample)
        expected = round(
            relative.coverage
            * relative.genome_sketch.genome_size
            / read_length
            * bases_explained
        )
        self.assertGreater(counts.sequence_abundance, 0)
        self.assertEqual(counts.sequence_abundance, float(expected))
        self.assertEqual(relative.taxonomic_abundance, counts.taxonomic_abundance)
        self.assertEqual(relative.coverage, counts.coverage)
