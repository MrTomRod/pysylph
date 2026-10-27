import io
import os
import pickle
import tempfile
import unittest
import importlib.resources
import pathlib

from pysylph import SampleSketch


def data_path(name):
    if hasattr(importlib.resources, "files"):
        return os.fspath(importlib.resources.files(__package__).joinpath(name))
    return os.fspath(pathlib.Path(__file__).parent / name)


class TestSampleSketchSerialization(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sketch = SampleSketch.load(data_path("k12.R1.sylsp"))
        cls.paired = SampleSketch.load(data_path("k12.paired.sylsp"))

    def assertSameSketch(self, actual, expected):
        self.assertIsInstance(actual, SampleSketch)
        self.assertEqual(actual.k, expected.k)
        self.assertEqual(actual.c, expected.c)
        self.assertEqual(actual.kmer_counts, expected.kmer_counts)

    def test_to_bytes_matches_dump(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "k12.R1.sylsp")
            self.sketch.dump(path)
            with open(path, "rb") as f:
                dumped = f.read()
        self.assertEqual(self.sketch.to_bytes(), dumped)

    def test_to_bytes_matches_file(self):
        with open(data_path("k12.R1.sylsp"), "rb") as f:
            self.assertEqual(self.sketch.to_bytes(), f.read())

    def test_from_bytes(self):
        sketch = SampleSketch.from_bytes(self.sketch.to_bytes())
        self.assertSameSketch(sketch, self.sketch)

    def test_from_bytes_other_samples(self):
        for expected in (self.sketch, self.paired):
            with self.subTest(expected=expected):
                sketch = SampleSketch.from_bytes(expected.to_bytes())
                self.assertSameSketch(sketch, expected)

    def test_roundtrip_buffer_types(self):
        data = self.sketch.to_bytes()
        for buffer in (bytes(data), bytearray(data), memoryview(data)):
            with self.subTest(type=type(buffer).__name__):
                sketch = SampleSketch.from_bytes(buffer)
                self.assertSameSketch(sketch, self.sketch)

    def test_bytes_io_interoperability(self):
        # `from_bytes` should read whatever `dump` writes...
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "k12.R1.sylsp")
            self.sketch.dump(path)
            with open(path, "rb") as f:
                sketch = SampleSketch.from_bytes(f.read())
        self.assertSameSketch(sketch, self.sketch)

        # ... and `load` should read whatever `to_bytes` produces
        sketch = SampleSketch.load(io.BytesIO(self.sketch.to_bytes()))
        self.assertSameSketch(sketch, self.sketch)

    def test_to_bytes_dump_to_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "copy.sylsp")
            with open(path, "wb") as f:
                f.write(self.sketch.to_bytes())
            self.assertSameSketch(SampleSketch.load(path), self.sketch)

    def test_from_bytes_invalid_data(self):
        for data in (b"", b"\x00", bytes(range(256)), b"not a sketch at all"):
            with self.subTest(data=data):
                self.assertRaises(ValueError, SampleSketch.from_bytes, data)

    def test_from_bytes_truncated_data(self):
        data = self.sketch.to_bytes()
        for size in (1, len(data) // 2, len(data) - 1):
            with self.subTest(size=size):
                self.assertRaises(ValueError, SampleSketch.from_bytes, data[:size])

    def test_from_bytes_invalid_type(self):
        self.assertRaises(TypeError, SampleSketch.from_bytes, 42)

    def test_pickle(self):
        sketch = pickle.loads(pickle.dumps(self.sketch))
        self.assertSameSketch(sketch, self.sketch)
