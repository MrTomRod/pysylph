import os
import unittest
import gzip
import random
import importlib.resources
from contextlib import nullcontext

from pysylph import Database, SampleSketch, Sketcher, Profiler


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
