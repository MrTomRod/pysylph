# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pysylph/compare/v0.1.3...HEAD

### Added
- `AniResult` properties `ani_ci`, `effective_lambda`, `lambda_status`,
  `lambda_ci`, `median_coverage`, `mean_coverage` and `containment`, exposing
  the remaining columns of the `sylph query` and `sylph profile` output.
- `min_contain` keyword argument to `Profiler` (default: 7), the minimum
  number of contained k-mers required to report a genome, as in `sylph` v1.0.0.
- `TwoStageDatabase` class to load, build (from a `Database`) and dump the
  two-stage databases (`.syl2db`) introduced in `sylph` v1.0.0.
- Support for `TwoStageDatabase` in `Profiler`, with a `screen_ani` keyword
  argument (default: 85) controlling the first-stage screen.
- `Profiler.database` is now a property returning either database type.

### Changed
- `AniResult` and `ProfileResult` now own their data, and genomes are matched
  to results by identity instead of by name, which was wrong when several
  genomes of a database share the same name.
- `Profiler.profile` and `Profiler.query` now release the GIL during the whole
  computation.
- `Database.load` and `SampleSketch.load` now release the GIL while decoding
  when given a path.
- Update `sylph` to `v1.0.0` (via a fork exposing its API, replacing the
  vendored copy of `sylph` internals).
- Default of `Profiler.min_number_kmers` lowered from 50 to 10 to match `sylph` v1.0.0.
- `Profiler.query` and `Profiler.profile` now raise `ValueError` when the sample
  was sketched with a different `k` than a database genome, or with a `c` greater
  than a database genome, instead of terminating the process.
- Read deduplication filter and bootstrap confidence intervals are now seeded
  as in `sylph` v1.0.0, making results reproducible.

### Fixed
- `Sketcher.sketch_single` never deduplicated reads because of an inverted
  read-length check.
- ANI values computed with a wrong value of `e` in the vendored
  `ani_from_lambda` (fixed in `sylph` v0.8.1); ANI values rise by about 0.06.
- `Profiler` treated `seq_id` as a fraction instead of a percentage, giving
  meaningless coverages when combined with `estimate_unknown`.


## [v0.1.3] - 2026-09-01
[v0.1.3]: https://github.com/althonos/pysylph/compare/v0.1.2...v0.1.3

### Added
- Wheels with Python 3.11 Limited API.

### Changed
- Update `sylph` to `v0.9.0`.
- Update `pyo3` to `v0.29.0`.


## [v0.1.2] - 2024-11-05
[v0.1.2]: https://github.com/althonos/pysylph/compare/v0.1.1...v0.1.2

### Added
- `Sketcher.sketch_paired` to sketch paired reads.
- `SampleSketch.kmer_counts` to access the k-mer counts of a sketch as a Python `dict`. 
- `deduplicated` and `fpr` parameters to `Sketcher` to control read deduplication.
- Package recipe on the Arch User Repository.


## [v0.1.1] - 2024-10-25
[v0.1.1]: https://github.com/althonos/pysylph/compare/v0.1.0...v0.1.1

### Added
- Type annotations to `pysylph` module.
- Sphinx documentation with API reference hosted on ReadTheDocs.

### Fixed
- Add missing `ProfileResult` class to the `pysylph.lib` module.

### Changed
- Implement reading sequence data using the buffer-protocol in `Sketcher` methods.


## [v0.1.0] - 2024-10-25
[v0.1.0]: https://github.com/althonos/pysylph/compare/3c6d23...v0.1.0

Initial release.
