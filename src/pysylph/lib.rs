extern crate bincode;
extern crate pyo3;
extern crate sylph;

use std::io::Read;
use std::os::unix::fs::FileExt;
use std::os::unix::fs::MetadataExt;
use std::sync::Arc;

use bincode::de::read::IoReader;
use pyo3::buffer::PyBuffer;
use pyo3::exceptions::PyIndexError;
use pyo3::exceptions::PyNotImplementedError;
use pyo3::exceptions::PyRuntimeError;
use pyo3::exceptions::PyTypeError;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedStr;
use pyo3::types::PyDict;
use pyo3::types::PyList;
use pyo3::types::PyString;
use pyo3::types::PyType;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use rayon::prelude::*;

mod pyfile;

// --- Base --------------------------------------------------------------------

/// A base sketch object.
#[pyclass(module = "pysylph.lib", frozen, subclass)]
pub struct Sketch;

#[pymethods]
impl Sketch {
    /// `int`: The subsampling rate this sketch was built with.
    #[getter]
    pub fn c(&self) -> PyResult<usize> {
        Err(PyNotImplementedError::new_err("not implemented"))
    }

    /// `int`: The k-mer size this sketch was built with.
    #[getter]
    pub fn k(&self) -> PyResult<usize> {
        Err(PyNotImplementedError::new_err("not implemented"))
    }
}

// --- GenomeSketch ------------------------------------------------------------

/// A (reference) genome sketch.
#[pyclass(module = "pysylph.lib", frozen, extends=Sketch)]
pub struct GenomeSketch {
    sketch: Arc<sylph::types::GenomeSketch>,
}

impl From<Arc<sylph::types::GenomeSketch>> for GenomeSketch {
    fn from(sketch: Arc<sylph::types::GenomeSketch>) -> Self {
        Self { sketch }
    }
}

impl From<sylph::types::GenomeSketch> for GenomeSketch {
    fn from(sketch: sylph::types::GenomeSketch) -> Self {
        Self::from(Arc::new(sketch))
    }
}

#[pymethods]
impl GenomeSketch {
    pub fn __repr__<'py>(slf: PyRef<'py, Self>) -> PyResult<String> {
        Ok(format!("<GenomeSketch name={:?}>", slf.sketch.file_name))
    }

    /// `str`: The name of the genome this sketch originates from.
    #[getter]
    pub fn name(&self) -> &str {
        &self.sketch.file_name
    }

    /// `str`: The description of the genome this sketch originates from.
    #[getter]
    pub fn description(&self) -> &str {
        &self.sketch.first_contig_name
    }

    /// `int`: The total number of nucleotides in the genome.
    #[getter]
    pub fn genome_size(&self) -> usize {
        self.sketch.gn_size
    }

    /// `int`: The subsampling rate this sketch was built with.
    #[getter]
    pub fn c(&self) -> usize {
        self.sketch.c
    }

    /// `int`: The k-mer size this sketch was built with.
    #[getter]
    pub fn k(&self) -> usize {
        self.sketch.k
    }

    /// `int`: The minimum spacing between k-mers used when sketching.
    #[getter]
    pub fn min_spacing(&self) -> usize {
        self.sketch.min_spacing
    }

    /// `list` of `int`: The list of k-mers extracted from the genome.
    #[getter]
    pub fn kmers<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(py, &self.sketch.genome_kmers)
    }
}

// --- Database ----------------------------------------------------------------

/// A sketch database.
///
/// Sylph databases are simply a sequence of `GenomeSketch` concatenated
/// together to be used for fast profiling or querying.
///
#[pyclass(module = "pysylph.lib", frozen)]
#[derive(Debug, Default)]
pub struct Database {
    sketches: Vec<Arc<sylph::types::GenomeSketch>>,
}

impl From<Vec<Arc<sylph::types::GenomeSketch>>> for Database {
    fn from(sketches: Vec<Arc<sylph::types::GenomeSketch>>) -> Self {
        Self { sketches }
    }
}

impl FromIterator<sylph::types::GenomeSketch> for Database {
    fn from_iter<T: IntoIterator<Item = sylph::types::GenomeSketch>>(iter: T) -> Self {
        let it = iter.into_iter();
        let sketches = it.map(Arc::new).collect::<Vec<_>>();
        Self::from(sketches)
    }
}

#[pymethods]
impl Database {
    #[new]
    #[pyo3(signature = (items = None))]
    pub fn __new__<'py>(items: Option<Bound<'py, PyAny>>) -> PyResult<Self> {
        let mut db = Self::default();
        if let Some(sketches) = items {
            for object in sketches.try_iter()? {
                let sketch: PyRef<'py, GenomeSketch> = object?.extract()?;
                db.sketches.push(sketch.sketch.clone());
            }
        }
        Ok(db)
    }

    pub fn __len__<'py>(slf: PyRef<'py, Self>) -> usize {
        slf.sketches.len()
    }

    pub fn __getitem__<'py>(slf: PyRef<'py, Self>, item: isize) -> PyResult<Py<GenomeSketch>> {
        let py = slf.py();
        let mut item_ = item;
        if item_ < 0 {
            item_ += slf.sketches.len() as isize;
        }

        if item_ < 0 || item_ >= slf.sketches.len() as isize {
            Err(PyIndexError::new_err(item))
        } else {
            Py::new(
                py,
                PyClassInitializer::from(Sketch)
                    .add_subclass(GenomeSketch::from(slf.sketches[item_ as usize].clone())),
            )
        }
    }

    /// Load a database from a file.
    ///
    /// Arguments:
    ///     file (`str`, `os.PathLike` or file-like handle): The path to the
    ///         file containing the database, or a file-like object open in
    ///         binary mode.
    ///
    /// Example:
    ///     >>> db = pysylph.Database.load("ecoli.syldb")
    ///     >>> db[0].name
    ///     'test_files/e.coli-K12.fasta.gz'
    ///
    #[classmethod]
    #[pyo3(signature = (file))]
    fn load<'py>(cls: &Bound<'py, PyType>, file: &Bound<'py, PyAny>) -> PyResult<Self> {
        let py = cls.py();
        // attempt to extract a path, or fall back to reader
        let result: Result<Vec<Arc<sylph::types::GenomeSketch>>, _> = if let Ok(s) = py
            .import(pyo3::intern!(py, "os"))?
            .call_method1(pyo3::intern!(py, "fsdecode"), (file,))
        {
            let path = s.cast::<PyString>()?;
            let f = std::fs::File::open(path.to_str()?)?;
            let reader = std::io::BufReader::new(f);
            bincode::deserialize_from(reader)
        } else {
            let f = self::pyfile::PyFileRead::from_ref(&file)?;
            let reader = std::io::BufReader::new(f);
            bincode::deserialize_from(reader)
        };
        // handle loader error
        match result {
            Ok(sketches) => Ok(Database::from(sketches)),
            Err(e) => match *e {
                bincode::ErrorKind::Io(io) => Err(io.into()),
                other => Err(PyValueError::new_err(format!(
                    "failed to load db: {:?}",
                    other
                ))),
            },
        }
    }

    /// Dump a database to a path.
    fn dump<'py>(slf: PyRef<'py, Self>, path: PyBackedStr) -> PyResult<()> {
        let f = std::fs::File::create(&*path)?;

        bincode::serialize_into(f, &slf.sketches).unwrap();
        Ok(())
    }
}

// --- Database Reader ---------------------------------------------------------

struct DatabaseReader<R: Read> {
    reader: bincode::de::Deserializer<
        IoReader<R>,
        bincode::config::WithOtherTrailing<
            bincode::config::WithOtherIntEncoding<
                bincode::config::DefaultOptions,
                bincode::config::FixintEncoding,
            >,
            bincode::config::AllowTrailing,
        >,
    >,
    length: usize,
}

impl<R: Read> DatabaseReader<R> {
    fn new(r: R) -> Result<Self, bincode::Error> {
        use bincode::Options;
        let mut reader = bincode::de::Deserializer::with_reader(
            r,
            bincode::config::DefaultOptions::new()
                .with_fixint_encoding()
                .allow_trailing_bytes(),
        );
        let length: usize = serde::Deserialize::deserialize(&mut reader).unwrap_or(0);
        Ok(Self { reader, length })
    }
}

impl<R: Read> std::iter::Iterator for DatabaseReader<R> {
    type Item = Result<sylph::types::GenomeSketch, bincode::Error>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.length > 0 {
            use serde::Deserialize;
            let item = sylph::types::GenomeSketch::deserialize(&mut self.reader).unwrap();
            self.length -= 1;
            Some(Ok(item))
        } else {
            None
        }
    }
}

/// A database file.
///
/// This class can be used to load `GenomeSketch` objects iteratively from a
/// database file instead of loading it entirely through the `Database.load`
/// class method.
///
/// Example:
///     >>> dbfile = pysylph.DatabaseFile("ecoli.syldb")
///     >>> for sketch in dbfile:
///     ...     print(sketch.name)
///     test_files/e.coli-K12.fasta.gz
///     test_files/e.coli-o157.fasta.gz
///     test_files/e.coli-EC590.fasta.gz
///
#[pyclass(module = "pysylph.lib")]
pub struct DatabaseFile {
    reader: DatabaseReader<std::io::BufReader<std::fs::File>>,
}

#[pymethods]
impl DatabaseFile {
    #[new]
    fn __new__(path: &str) -> PyResult<Self> {
        let f = std::fs::File::open(path).map(std::io::BufReader::new)?;
        let reader = DatabaseReader::new(f).unwrap();
        Ok(Self { reader })
    }

    fn __len__(&self) -> usize {
        self.reader.length
    }

    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }

    fn __next__<'py>(mut slf: PyRefMut<'py, Self>) -> PyResult<Option<Py<GenomeSketch>>> {
        let py = slf.py();
        match slf.reader.next() {
            Some(Err(e)) => Err(PyRuntimeError::new_err(e.to_string())),
            None => Ok(None),
            Some(Ok(item)) => Some(Py::new(
                py,
                PyClassInitializer::from(Sketch).add_subclass(GenomeSketch::from(item)),
            ))
            .transpose(),
        }
    }

    // fn read_all(&mut self) -> Database {
    //     let mut db = Database::from(Vec::with_capacity(self.reader.length));
    //     for sketch in &mut self.reader {
    //         db.sketches.push(Arc::new(sketch.unwrap()));
    //     }
    //     let vec = Vec::<Arc<sylph::types::GenomeSketch>>::deserialize(&mut self.reader.reader).unwrap();
    //     Database::from(vec)
    // }

    // fn read(&mut self) -> GenomeSketch {
    //     use serde::de::Deserialize;
    //     let l: usize = serde::Deserialize::deserialize(&mut d).unwrap();
    //     // println!("{:?}", d);
    //     unimplemented!()
    // }
}

// --- TwoStageDatabase --------------------------------------------------------

/// Convert an I/O error raised while loading a database to a Python error.
fn twostage_error(e: std::io::Error) -> PyErr {
    use std::io::ErrorKind;
    match e.kind() {
        ErrorKind::InvalidData | ErrorKind::UnexpectedEof | ErrorKind::InvalidInput => {
            PyValueError::new_err(format!("invalid two-stage database: {}", e))
        }
        _ => e.into(),
    }
}

/// Run some ``sylph`` code, turning panics caused by malformed data into errors.
fn guard_panic<T>(f: impl FnOnce() -> std::io::Result<T>) -> std::io::Result<T> {
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)).unwrap_or_else(|_| {
        Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "malformed two-stage database",
        ))
    })
}

/// A two-stage sketch database.
///
/// Two-stage databases (``.syl2db``) were introduced in ``sylph`` v1.0.0.
/// They store a sparse index of every genome, used to quickly screen
/// samples, along with the dense k-mers of each genome, which are only
/// decoded for the genomes that pass the screen. This makes profiling
/// large databases of large genomes much faster, with little memory.
///
/// Use `TwoStageDatabase.from_database` to build one from a `Database`, or
/// `TwoStageDatabase.load` to load one built with ``sylph``. Two-stage
/// databases can be given to `Profiler` in place of a `Database`.
///
/// .. versionadded:: 0.2.0
///
#[pyclass(module = "pysylph.lib", frozen)]
pub struct TwoStageDatabase {
    db: sylph::twostage_db::TwoStageDb,
    /// For databases read lazily from a file: the file, kept open so that it
    /// can still be copied if the path is removed or replaced. In-memory
    /// databases have no handle, and are dumped from `TwoStageDb::as_bytes`.
    file: Option<std::fs::File>,
}

#[pymethods]
impl TwoStageDatabase {
    /// `int`: The subsampling rate of the dense k-mers of each genome.
    #[getter]
    pub fn c(&self) -> usize {
        self.db.c
    }

    /// `int`: The k-mer size the genomes were sketched with.
    #[getter]
    pub fn k(&self) -> usize {
        self.db.k
    }

    /// `int`: The subsampling rate of the screening index.
    ///
    /// This can be lower than the rate requested when building the database
    /// if small genomes needed a denser screen.
    #[getter]
    pub fn screen_c(&self) -> usize {
        self.db.screen_c
    }

    pub fn __len__(&self) -> usize {
        self.db.len()
    }

    pub fn __getitem__<'py>(slf: PyRef<'py, Self>, item: isize) -> PyResult<Py<GenomeSketch>> {
        let py = slf.py();
        let len = slf.db.len() as isize;
        let index = if item < 0 { item + len } else { item };
        if index < 0 || index >= len {
            return Err(PyIndexError::new_err(item));
        }
        let this: &Self = &slf;
        let sketch = py
            .detach(|| guard_panic(|| this.db.decode_dense(index as u32)))
            .map_err(twostage_error)?;
        Py::new(
            py,
            PyClassInitializer::from(Sketch).add_subclass(GenomeSketch::from(sketch)),
        )
    }

    /// Load a two-stage database.
    ///
    /// Arguments:
    ///     file (`str`, `os.PathLike` or file-like handle): The path to the
    ///         file containing the database, or a file-like object open in
    ///         binary mode.
    ///
    /// When given a path, only the screening index is loaded in memory, and
    /// the dense k-mers of genomes are read from the file when needed. When
    /// given a file-like object, the whole database is loaded in memory.
    ///
    /// Raises:
    ///     `ValueError`: When the file is not a valid two-stage database.
    ///
    #[classmethod]
    #[pyo3(signature = (file))]
    fn load<'py>(cls: &Bound<'py, PyType>, file: &Bound<'py, PyAny>) -> PyResult<Self> {
        let py = cls.py();
        // attempt to extract a path, or fall back to reader
        if let Ok(s) = py
            .import(pyo3::intern!(py, "os"))?
            .call_method1(pyo3::intern!(py, "fsdecode"), (file,))
        {
            let path = s.cast::<PyString>()?.to_str()?.to_string();
            let handle = std::fs::File::open(&path)?;
            let db = py
                .detach(|| {
                    guard_panic(|| sylph::twostage_db::open_file(&path))
                })
                .map_err(twostage_error)?;
            Ok(Self {
                db,
                file: Some(handle),
            })
        } else {
            let mut reader = self::pyfile::PyFileRead::from_ref(&file)?;
            let mut bytes = Vec::new();
            reader.read_to_end(&mut bytes)?;
            let db = py
                .detach(move || guard_panic(move || sylph::twostage_db::from_vec(bytes)))
                .map_err(twostage_error)?;
            Ok(Self { db, file: None })
        }
    }

    /// Build a two-stage database from a database of genome sketches.
    ///
    /// Arguments:
    ///     database (`~pysylph.Database`): The database to convert. All of
    ///         its sketches must have the same ``c`` and ``k``, and have been
    ///         sketched with ``profiling=True``.
    ///     screen_c (`int`): The subsampling rate of the screening index,
    ///         which must be greater or equal to the ``c`` of the genomes.
    ///     min_sparse_kmers (`int`): The minimum number of k-mers of each
    ///         genome in the screening index. Genomes with less k-mers at the
    ///         ``screen_c`` rate use a denser, genome-specific rate.
    ///     min_contain (`int`): Genomes with fewer dense k-mers than this
    ///         are dropped from the database, since they could not be
    ///         reported by a `Profiler` with the same ``min_contain``. As a
    ///         result, the database can have less genomes than ``database``.
    ///
    /// The database is built in memory.
    ///
    /// Raises:
    ///     `ValueError`: When the database is empty, has inconsistent
    ///         ``c`` or ``k``, has genomes sketched without profiling
    ///         k-mers, or when a parameter is out of range.
    ///
    #[classmethod]
    #[pyo3(signature = (database, *, screen_c = sylph::constants::SCREEN_C_DEFAULT, min_sparse_kmers = sylph::constants::SPARSE_TARGET_MIN_DEFAULT, min_contain = 7))]
    fn from_database<'py>(
        cls: &Bound<'py, PyType>,
        database: &Database,
        screen_c: usize,
        min_sparse_kmers: usize,
        min_contain: usize,
    ) -> PyResult<Self> {
        let py = cls.py();
        let sketches = &database.sketches;
        let first = match sketches.first() {
            Some(first) => first,
            None => return Err(PyValueError::new_err("cannot convert an empty database")),
        };
        for sketch in sketches.iter() {
            if sketch.c != first.c || sketch.k != first.k {
                return Err(PyValueError::new_err(format!(
                    "inconsistent c or k in database: genome {} has c={} and k={}, but genome {} has c={} and k={}",
                    first.file_name, first.c, first.k, sketch.file_name, sketch.c, sketch.k,
                )));
            }
            if sketch.pseudotax_tracked_nonused_kmers.is_none() {
                return Err(PyValueError::new_err(format!(
                    "genome {} was sketched without profiling k-mers, which two-stage databases require",
                    sketch.file_name
                )));
            }
        }
        if screen_c < first.c {
            return Err(PyValueError::new_err(format!(
                "screen_c ({}) must be greater or equal to the database c ({})",
                screen_c, first.c
            )));
        }
        if min_sparse_kmers < 1 {
            return Err(PyValueError::new_err("min_sparse_kmers must be at least 1"));
        }

        let result = py.detach(|| {
            guard_panic(|| {
                let mut bytes = Vec::new();
                sylph::twostage_db::write_two_stage_db(
                    &mut bytes,
                    &sketches[..],
                    screen_c,
                    min_sparse_kmers,
                    min_contain,
                )?;
                sylph::twostage_db::from_vec(bytes)
            })
        });
        match result {
            Ok(db) => Ok(Self { db, file: None }),
            Err(e) => Err(PyValueError::new_err(format!(
                "failed to build two-stage database: {}",
                e
            ))),
        }
    }

    /// Dump a two-stage database to a path.
    fn dump<'py>(slf: PyRef<'py, Self>, path: &Bound<'py, PyAny>) -> PyResult<()> {
        let py = slf.py();
        let path = py
            .import(pyo3::intern!(py, "os"))?
            .call_method1(pyo3::intern!(py, "fsdecode"), (path,))?
            .cast::<PyString>()?
            .to_str()?
            .to_string();
        let this: &Self = &slf;
        py.detach(|| match (this.db.as_bytes(), &this.file) {
            (Some(bytes), _) => std::fs::write(&path, bytes),
            (None, Some(src)) => {
                // do not truncate the database if dumping it over itself
                let meta = src.metadata()?;
                if let Ok(dst) = std::fs::metadata(&path) {
                    if dst.dev() == meta.dev() && dst.ino() == meta.ino() {
                        return Ok(());
                    }
                }
                let mut dst = std::fs::File::create(&path)?;
                let mut buffer = vec![0u8; 1 << 20];
                let mut offset = 0u64;
                while offset < meta.len() {
                    let n = usize::min(buffer.len(), (meta.len() - offset) as usize);
                    src.read_exact_at(&mut buffer[..n], offset)?;
                    std::io::Write::write_all(&mut dst, &buffer[..n])?;
                    offset += n as u64;
                }
                Ok(())
            }
            (None, None) => unreachable!("file-backed database without a file handle"),
        })?;
        Ok(())
    }
}

// --- SampleSketch ----------------------------------------------------------

/// A (query) sample sketch.
#[pyclass(module = "pysylph.lib", frozen, extends=Sketch)]
pub struct SampleSketch {
    sketch: sylph::types::SequencesSketch,
}

impl From<sylph::types::SequencesSketch> for SampleSketch {
    fn from(sketch: sylph::types::SequencesSketch) -> Self {
        Self { sketch }
    }
}

#[pymethods]
impl SampleSketch {
    pub fn __repr__<'py>(slf: PyRef<'py, Self>) -> PyResult<String> {
        Ok(format!("<SampleSketch name={:?}>", slf.sketch.file_name))
    }

    /// `int`: The subsampling rate this sketch was built with.
    #[getter]
    pub fn c(&self) -> usize {
        self.sketch.c
    }

    /// `int`: The k-mer size this sketch was built with.
    #[getter]
    pub fn k(&self) -> usize {
        self.sketch.k
    }

    /// `dict`: A counter of k-mers in this sketch.
    #[getter]
    pub fn kmer_counts<'py>(slf: PyRef<'py, Self>) -> PyResult<Bound<'py, PyDict>> {
        let py = slf.py();
        let mapping = PyDict::new(py);
        for (kmer, count) in slf.sketch.kmer_counts.iter() {
            mapping.set_item(kmer, count)?;
        }
        Ok(mapping)
    }

    /// Load a sequence sketch from a path.
    #[classmethod]
    #[pyo3(signature = (file))]
    fn load<'py>(cls: &Bound<'py, PyType>, file: &Bound<'py, PyAny>) -> PyResult<Py<Self>> {
        let py = cls.py();
        // attempt to extract a path, or fall back to reader
        let result: Result<sylph::types::SequencesSketch, _> = if let Ok(s) = py
            .import(pyo3::intern!(py, "os"))?
            .call_method1(pyo3::intern!(py, "fsdecode"), (file,))
        {
            let path = s.cast::<PyString>()?;
            let f = std::fs::File::open(path.to_str()?)?;
            let reader = std::io::BufReader::new(f);
            bincode::deserialize_from(reader)
        } else {
            let f = self::pyfile::PyFileRead::from_ref(&file)?;
            let reader = std::io::BufReader::new(f);
            bincode::deserialize_from(reader)
        };
        // handle error
        match result {
            Ok(sketch) => Py::new(
                py,
                PyClassInitializer::from(Sketch).add_subclass(Self::from(sketch)),
            ),
            Err(e) => match *e {
                bincode::ErrorKind::Io(io) => Err(io.into()),
                other => Err(PyValueError::new_err(format!(
                    "failed to load db: {:?}",
                    other
                ))),
            },
        }
    }

    /// Dump a sequence sketch to a path.
    fn dump<'py>(slf: PyRef<'py, Self>, path: PyBackedStr) -> PyResult<()> {
        let f = std::fs::File::create(&*path)?;

        bincode::serialize_into(f, &slf.sketch).unwrap();
        Ok(())
    }
}

// --- ANIResult ---------------------------------------------------------------

/// An querying result.
#[pyclass(module = "pysylph.lib", frozen, subclass)]
pub struct AniResult {
    final_est_ani: f64,
    naive_ani: f64,
    final_est_cov: f64,
    mean_cov: f64,
    median_cov: f64,
    containment_index: (usize, usize),
    lambda: sylph::types::AdjustStatus,
    ani_ci: (Option<f64>, Option<f64>),
    lambda_ci: (Option<f64>, Option<f64>),
    genome: Py<GenomeSketch>,
}

#[pymethods]
impl AniResult {
    pub fn __repr__<'py>(slf: PyRef<'py, Self>) -> PyResult<String> {
        let py = slf.py();
        let genome = slf.genome.bind(py).borrow();
        Ok(format!(
            "<AniResult genome={:?} ani={:?}>",
            genome.sketch.file_name, slf.final_est_ani
        ))
    }

    /// `~pysylph.GenomeSketch`: A reference to the genome sketch.
    #[getter]
    fn genome_sketch<'py>(slf: PyRef<'py, Self>) -> Py<GenomeSketch> {
        let py = slf.py();
        slf.genome.clone_ref(py)
    }

    /// `float`: The coverage-corrected containment ANI, as a percentage.
    #[getter]
    fn ani<'py>(slf: PyRef<'py, Self>) -> f64 {
        f64::min(slf.final_est_ani * 100.0, 100.0)
    }

    /// `float`: The containment ANI without adjustment, as a percentage.
    #[getter]
    fn ani_naive<'py>(slf: PyRef<'py, Self>) -> f64 {
        slf.naive_ani * 100.0
    }

    /// `float`: An estimate of the effective, or the true coverage.
    ///
    /// This is the ``Eff_cov`` column of ``sylph``, or the ``True_cov`` column
    /// when the `Profiler` was created with ``estimate_unknown=True``.
    #[getter]
    fn coverage<'py>(slf: PyRef<'py, Self>) -> f64 {
        slf.final_est_cov
    }

    /// `tuple` of `float`, or `None`: The 5-95 percentile confidence interval
    /// of the adjusted ANI, as percentages.
    ///
    /// This is `None` when the interval could not be computed (``NA-NA`` in
    /// the ``sylph`` output).
    #[getter]
    fn ani_ci<'py>(slf: PyRef<'py, Self>) -> Option<(f64, f64)> {
        match slf.ani_ci {
            (Some(low), Some(high)) => Some((low * 100.0, high * 100.0)),
            _ => None,
        }
    }

    /// `float`, or `None`: The effective k-mer coverage lambda estimated by
    /// the Poisson model (``Eff_lambda``).
    ///
    /// This is `None` when the estimate is out of the supported range; see
    /// `AniResult.lambda_status` for the reason.
    #[getter]
    fn effective_lambda<'py>(slf: PyRef<'py, Self>) -> Option<f64> {
        match slf.lambda {
            sylph::types::AdjustStatus::Lambda(x) => Some(x),
            _ => None,
        }
    }

    /// `str`: The status of the lambda estimate.
    ///
    /// One of ``"OK"`` when `AniResult.effective_lambda` is available, or
    /// ``"HIGH"`` or ``"LOW"`` when the coverage was too high or too low to
    /// be estimated (as printed in the ``Eff_lambda`` column of ``sylph``).
    #[getter]
    fn lambda_status<'py>(slf: PyRef<'py, Self>) -> &'static str {
        match slf.lambda {
            sylph::types::AdjustStatus::Lambda(_) => "OK",
            sylph::types::AdjustStatus::High => "HIGH",
            sylph::types::AdjustStatus::Low => "LOW",
        }
    }

    /// `tuple` of `float`, or `None`: The 5-95 percentile confidence interval
    /// of the effective lambda (``Lambda_5-95_percentile``).
    ///
    /// This is `None` when the interval could not be computed.
    #[getter]
    fn lambda_ci<'py>(slf: PyRef<'py, Self>) -> Option<(f64, f64)> {
        match slf.lambda_ci {
            (Some(low), Some(high)) => Some((low, high)),
            _ => None,
        }
    }

    /// `float`: The median k-mer coverage of the genome (``Median_cov``).
    #[getter]
    fn median_coverage<'py>(slf: PyRef<'py, Self>) -> f64 {
        slf.median_cov
    }

    /// `float`: The mean k-mer coverage over the covered k-mers, i.e. those
    /// with a coverage of at least one (``Mean_cov_geq1``).
    #[getter]
    fn mean_coverage<'py>(slf: PyRef<'py, Self>) -> f64 {
        slf.mean_cov
    }

    /// `tuple` of `int`: The number of genome k-mers found in the sample, and
    /// the total number of genome k-mers (``Containment_ind``).
    #[getter]
    fn containment<'py>(slf: PyRef<'py, Self>) -> (usize, usize) {
        slf.containment_index
    }
}

/// A profiling result.
#[pyclass(module = "pysylph.lib", frozen, extends=AniResult)]
pub struct ProfileResult {
    sequence_abundance: f64,
    taxonomic_abundance: f64,
    kmers_reassigned: usize,
}

#[pymethods]
impl ProfileResult {
    pub fn __repr__<'py>(slf: PyRef<'py, Self>) -> PyResult<String> {
        let py = slf.py();
        let genome = slf.as_super().genome.bind(py).borrow();
        Ok(format!(
            "<ProfileResult genome={:?} ani={:?}>",
            genome.sketch.file_name,
            slf.as_super().final_est_ani
        ))
    }

    /// `float`: The normalized sequence abundance, as a percentage.
    ///
    /// This is the *percentage of reads* assigned to each genome, similar
    /// to how Kraken computes abundance.
    #[getter]
    fn sequence_abundance<'py>(slf: PyRef<'py, Self>) -> f64 {
        slf.sequence_abundance
    }

    /// `float`: The normalized taxonomic abundance, as a percentage. Coverage-normalized - same as MetaPhlAn abundance
    ///
    /// This is the *coverage-normalized* abundanced, similar to how
    /// MetaPhlAn computes abundance.
    #[getter]
    fn taxonomic_abundance<'py>(slf: PyRef<'py, Self>) -> f64 {
        slf.taxonomic_abundance
    }

    /// `int`: The number of k-mers reassigned away from the genome.
    #[getter]
    fn kmers_reassigned<'py>(slf: PyRef<'py, Self>) -> usize {
        slf.kmers_reassigned
    }
}

// --- Sketcher ----------------------------------------------------------------

enum SequenceData {
    #[allow(dead_code)]
    Buffer(PyBuffer<u8>, &'static [u8]),
    BackedStr(PyBackedStr),
}

impl SequenceData {
    fn as_bytes<'py>(&self) -> &[u8] {
        match self {
            SequenceData::BackedStr(s) => s.as_bytes(),
            SequenceData::Buffer(_, s) => s,
        }
    }
}

impl<'py, 'a> FromPyObject<'py, 'a> for SequenceData {
    type Error = PyErr;
    fn extract(obj: Borrowed<'py, 'a, PyAny>) -> PyResult<Self> {
        let py = obj.py();
        if let Ok(s) = obj.extract::<PyBackedStr>() {
            Ok(SequenceData::BackedStr(s))
        } else {
            match PyBuffer::get(&obj) {
                Ok(buffer) => {
                    if !buffer.is_c_contiguous() {
                        Err(PyValueError::new_err("expected C-contiguous buffer"))
                    } else if buffer.shape().len() != 1 {
                        Err(PyValueError::new_err("expected buffer of dimension 1"))
                    } else {
                        // NOTE(@althonos): Yes, this is unsafe because the Python
                        //                  object may be modified in parallel. We're
                        //                  gonna have to trust the user doesn't do
                        //                  weird shit with the sequence, but that's
                        //                  a trade-off so we can avoid copy data
                        //                  from the Python heap (or other Python
                        //                  extension objects).
                        let x = buffer
                            .as_slice(obj.py())
                            .ok_or_else(|| PyValueError::new_err("invalid buffer"))?;
                        let s =
                            unsafe { std::slice::from_raw_parts(x.as_ptr() as *const u8, x.len()) };
                        Ok(SequenceData::Buffer(buffer, s))
                    }
                }
                Err(t) => {
                    let e = PyTypeError::new_err("expected string or byte buffer");
                    e.set_cause(py, Some(t));
                    Err(e)
                }
            }
        }
    }
}

/// A ``sylph`` sketcher.
#[pyclass(module = "pysylph.lib", frozen)]
pub struct Sketcher {
    c: usize,
    k: usize,
    min_spacing: usize,
    deduplicate: bool,
    fpr: f64,
}

#[pymethods]
impl Sketcher {
    #[new]
    #[pyo3(signature = (*, c = 200, k = 31, deduplicate = true, fpr = sylph::constants::DEFAULT_FPR))]
    pub fn __new__(c: usize, k: usize, deduplicate: bool, fpr: f64) -> PyResult<Self> {
        if k != 21 && k != 31 {
            return Err(PyValueError::new_err(format!(
                "invalid k: expected 21 or 31, got {}",
                k
            )));
        }
        Ok(Self {
            c,
            k,
            min_spacing: 30,
            deduplicate,
            fpr,
        })
    }

    /// Sketch a genome.
    ///
    /// Arguments:
    ///     name (`str`): The name of the genome to sketch. In the ``sylph``
    ///         binary, this is set as the filename of the genome file.
    ///     contigs (iterable of `str` or byte buffer): The genome contigs
    ///         to be sketched.
    ///     profiling (`bool`): Set to `False` to disable the tracking of
    ///         unsued k-mers, which is required for profiling. This will
    ///         prevent the sketch from being used with `Profiler.profile`.
    ///
    /// Returns:
    ///     `~pysylph.GenomeSketch`: The sketched genome.
    ///
    #[pyo3(signature = (name, contigs, profiling=true))]
    fn sketch_genome<'py>(
        slf: PyRef<'py, Self>,
        name: String,
        contigs: Bound<'py, PyAny>,
        profiling: bool,
    ) -> PyResult<Py<GenomeSketch>> {
        let py = slf.py();

        let mut gsketch = sylph::types::GenomeSketch::default();
        gsketch.min_spacing = slf.min_spacing;
        gsketch.c = slf.c;
        gsketch.k = slf.k;
        gsketch.file_name = name;
        if profiling {
            gsketch.pseudotax_tracked_nonused_kmers = Some(Vec::new());
        }

        // extract records
        let sequences = contigs
            .try_iter()?
            .map(|r| r?.extract::<SequenceData>())
            .collect::<PyResult<Vec<_>>>()?;

        // sketch all records while allowing parallel code
        py.detach(|| {
            // extract candidate kmers
            let mut markers = Vec::new();
            for (index, sequence) in sequences.iter().enumerate() {
                sylph::sketch::extract_markers_positions(
                    sequence.as_bytes(),
                    &mut markers,
                    gsketch.c,
                    gsketch.k,
                    index,
                );
                gsketch.gn_size += sequence.as_bytes().len();
            }

            // split duplicate / unique kmers
            let mut kmer_set = sylph::types::MMHashSet::default();
            let mut duplicate_set = sylph::types::MMHashSet::default();
            markers.sort(); // NB(@althonos): is this necessary here?
            for (_, _, km) in markers.iter() {
                if !kmer_set.insert(km) {
                    duplicate_set.insert(km);
                }
            }

            // record kmers
            let mut last_pos = 0;
            let mut last_contig = 0;
            for &(contig, pos, km) in markers.iter() {
                if !duplicate_set.contains(&km) {
                    if last_pos == 0
                        || last_contig != contig
                        || pos > gsketch.min_spacing + last_pos
                    {
                        gsketch.genome_kmers.push(km);
                        last_contig = contig;
                        last_pos = pos;
                    } else if let Some(kmers) = &mut gsketch.pseudotax_tracked_nonused_kmers {
                        kmers.push(km);
                    }
                }
            }
        });

        Py::new(
            py,
            PyClassInitializer::from(Sketch).add_subclass(GenomeSketch::from(gsketch)),
        )
    }

    /// Sketch a sample composed of single-read sequences.
    ///
    /// Arguments:
    ///     name (`str`): The name of the sample to sketch. In the ``sylph``
    ///         binary, this is set as the filename of the genome file.
    ///     reads (iterable of `str` or byte buffer): The sample reads
    ///         to be sketched.
    ///
    /// Returns:
    ///     `~pysylph.SampleSketch`: The sketched sample.
    ///
    #[pyo3(signature = (name, reads))]
    fn sketch_single<'py>(
        slf: PyRef<'py, Self>,
        name: String,
        reads: Bound<'py, PyAny>,
    ) -> PyResult<Py<SampleSketch>> {
        let py = slf.py();

        let mut kmer_map = std::collections::HashMap::default();
        let mut mean_read_length = 0.;
        let mut counter = 0usize;
        let mut kmer_to_pair_table = fxhash::FxHashSet::default();
        let mut num_dup_removed = 0;

        for result in reads.try_iter()? {
            let read = result?.extract::<SequenceData>()?;
            let seq = read.as_bytes();

            let mut vec = vec![];
            let kmer_pair = if seq.len() > 400 {
                None
            } else {
                sylph::sketch::pair_kmer_single(seq)
            };
            sylph::sketch::extract_markers(&seq, &mut vec, slf.c, slf.k);
            for km in vec {
                sylph::sketch::dup_removal_lsh_full_exact(
                    &mut kmer_map,
                    &mut kmer_to_pair_table,
                    &km,
                    kmer_pair,
                    &mut num_dup_removed,
                    !slf.deduplicate,
                    Some(sylph::constants::MAX_DEDUP_COUNT),
                );
            }
            //moving average
            counter += 1;
            mean_read_length += ((seq.len() as f64) - mean_read_length) / counter as f64;
        }

        let sketch = sylph::types::SequencesSketch {
            kmer_counts: kmer_map,
            file_name: name,
            c: slf.c,
            k: slf.k,
            paired: false,
            sample_name: None,
            mean_read_length,
        };

        Py::new(
            py,
            PyClassInitializer::new(SampleSketch::from(sketch), PyClassInitializer::from(Sketch)),
        )
    }

    /// Sketch a sample composed of paired reads.
    ///
    /// Arguments:
    ///     name (`str`): The name of the sample to sketch. In the ``sylph``
    ///         binary, this is set as the filename of the genome file.
    ///     r1 (iterable of `str` or byte buffer): The first sequences
    ///         for paired end reads to be sketched.
    ///     r2 (iterable of `str` or byte buffer): The second sequences
    ///         for paired end reads to be sketched.
    ///
    /// Returns:
    ///     `~pysylph.SampleSketch`: The sketched sample.
    ///
    /// .. versionadded:: 0.1.2
    ///
    #[pyo3(signature = (name, r1, r2))]
    fn sketch_paired<'py>(
        slf: PyRef<'py, Self>,
        name: String,
        r1: Bound<'py, PyAny>,
        r2: Bound<'py, PyAny>,
    ) -> PyResult<Py<SampleSketch>> {
        let py = slf.py();

        let mut kmer_map = std::collections::HashMap::default();
        let mut mean_read_length = 0.;
        let mut counter = 0usize;
        let mut kmer_pair_set = fxhash::FxHashSet::default();
        let mut num_dup_removed = 0;

        let mut kmer_pair_set_approx = match slf.fpr {
            0.0 => None,
            x => Some(
                scalable_cuckoo_filter::ScalableCuckooFilterBuilder::new()
                    .initial_capacity(1_000_000_0)
                    .false_positive_probability(x)
                    .hasher(fxhash::FxHasher::default())
                    .rng(SmallRng::seed_from_u64(sylph::constants::DEFAULT_RNG_SEED))
                    .finish(),
            ),
        };

        for (res1, res2) in r1.try_iter()?.zip(r2.try_iter()?) {
            let r1 = res1?.extract::<SequenceData>()?;
            let r2 = res2?.extract::<SequenceData>()?;

            let seq1 = r1.as_bytes();
            let seq2 = r2.as_bytes();

            let mut vec1 = vec![];
            let mut vec2 = vec![];

            sylph::sketch::extract_markers(&seq1, &mut vec1, slf.c, slf.k);
            sylph::sketch::extract_markers(&seq2, &mut vec2, slf.c, slf.k);
            let kmer_pair = sylph::sketch::pair_kmer(&seq1, &seq2);

            //moving average
            counter += 1;
            mean_read_length += ((seq1.len() as f64) - mean_read_length) / counter as f64;

            for km in vec1.iter() {
                if let Some(mut approx) = kmer_pair_set_approx.as_mut() {
                    sylph::sketch::dup_removal_lsh_full(
                        &mut kmer_map,
                        &mut approx,
                        km,
                        kmer_pair,
                        &mut num_dup_removed,
                        !slf.deduplicate,
                    );
                } else {
                    sylph::sketch::dup_removal_lsh_full_exact(
                        &mut kmer_map,
                        &mut kmer_pair_set,
                        km,
                        kmer_pair,
                        &mut num_dup_removed,
                        !slf.deduplicate,
                        None,
                    );
                }
            }
            for km in vec2.iter() {
                if vec1.contains(km) {
                    continue;
                }
                if let Some(mut approx) = kmer_pair_set_approx.as_mut() {
                    sylph::sketch::dup_removal_lsh_full(
                        &mut kmer_map,
                        &mut approx,
                        km,
                        kmer_pair,
                        &mut num_dup_removed,
                        !slf.deduplicate,
                    );
                } else {
                    sylph::sketch::dup_removal_lsh_full_exact(
                        &mut kmer_map,
                        &mut kmer_pair_set,
                        km,
                        kmer_pair,
                        &mut num_dup_removed,
                        !slf.deduplicate,
                        None,
                    );
                }
            }
        }

        let sketch = sylph::types::SequencesSketch {
            kmer_counts: kmer_map,
            file_name: name,
            c: slf.c,
            k: slf.k,
            paired: true,
            sample_name: None,
            mean_read_length,
        };
        Py::new(
            py,
            PyClassInitializer::new(SampleSketch::from(sketch), PyClassInitializer::from(Sketch)),
        )
    }
}

// --- Functions ---------------------------------------------------------------

/// Check that the database and sample sketches can be compared.
///
/// ``sylph`` terminates the process when this does not hold, so make
/// sure to raise a Python exception before reaching that point.
fn check_compatible(database: &ProfilerDatabase, sample: &SampleSketch) -> PyResult<()> {
    match database {
        ProfilerDatabase::Plain(database) => {
            for genome in database.get().sketches.iter() {
                if genome.k != sample.sketch.k {
                    return Err(PyValueError::new_err(format!(
                        "k parameter for sample ({}) != k parameter for genome {} ({})",
                        sample.sketch.k, genome.file_name, genome.k,
                    )));
                }
                if genome.c < sample.sketch.c {
                    return Err(PyValueError::new_err(format!(
                        "c parameter for sample ({}) > c parameter for genome {} ({})",
                        sample.sketch.c, genome.file_name, genome.c,
                    )));
                }
            }
        }
        ProfilerDatabase::TwoStage(database) => {
            let db = &database.get().db;
            if db.k != sample.sketch.k {
                return Err(PyValueError::new_err(format!(
                    "k parameter for sample ({}) != k parameter for database ({})",
                    sample.sketch.k, db.k,
                )));
            }
            if db.c < sample.sketch.c {
                return Err(PyValueError::new_err(format!(
                    "c parameter for sample ({}) > c parameter for database ({})",
                    sample.sketch.c, db.c,
                )));
            }
        }
    }
    Ok(())
}

/// A hit of a sample against a genome, owning all the data needed to expose it.
struct Hit {
    genome: Arc<sylph::types::GenomeSketch>,
    final_est_ani: f64,
    naive_ani: f64,
    final_est_cov: f64,
    mean_cov: f64,
    median_cov: f64,
    containment_index: (usize, usize),
    lambda: sylph::types::AdjustStatus,
    ani_ci: (Option<f64>, Option<f64>),
    lambda_ci: (Option<f64>, Option<f64>),
    rel_abund: Option<f64>,
    seq_abund: Option<f64>,
    kmers_lost: Option<usize>,
}

impl Hit {
    fn new(
        result: &sylph::types::AniResult,
        genomes: &std::collections::HashMap<*const sylph::types::GenomeSketch, &Arc<sylph::types::GenomeSketch>>,
    ) -> Self {
        let genome = genomes[&(result.genome_sketch as *const sylph::types::GenomeSketch)];
        Self {
            genome: genome.clone(),
            final_est_ani: result.final_est_ani,
            naive_ani: result.naive_ani,
            final_est_cov: result.final_est_cov,
            mean_cov: result.mean_cov,
            median_cov: result.median_cov,
            containment_index: result.containment_index,
            lambda: match result.lambda {
                sylph::types::AdjustStatus::Lambda(x) => sylph::types::AdjustStatus::Lambda(x),
                sylph::types::AdjustStatus::High => sylph::types::AdjustStatus::High,
                sylph::types::AdjustStatus::Low => sylph::types::AdjustStatus::Low,
            },
            ani_ci: result.ani_ci,
            lambda_ci: result.lambda_ci,
            rel_abund: result.rel_abund,
            seq_abund: result.seq_abund,
            kmers_lost: result.kmers_lost,
        }
    }

    fn into_ani_result(self, py: Python<'_>) -> PyResult<Py<AniResult>> {
        let genome = Py::new(
            py,
            PyClassInitializer::from(Sketch).add_subclass(GenomeSketch::from(self.genome)),
        )?;
        Py::new(
            py,
            PyClassInitializer::from(AniResult {
                final_est_ani: self.final_est_ani,
                naive_ani: self.naive_ani,
                final_est_cov: self.final_est_cov,
                mean_cov: self.mean_cov,
                median_cov: self.median_cov,
                containment_index: self.containment_index,
                lambda: self.lambda,
                ani_ci: self.ani_ci,
                lambda_ci: self.lambda_ci,
                genome,
            }),
        )
    }

    fn into_profile_result(self, py: Python<'_>) -> PyResult<Py<ProfileResult>> {
        let profile = ProfileResult {
            sequence_abundance: self.seq_abund.ok_or_else(|| {
                PyRuntimeError::new_err("profile result without sequence abundance")
            })?,
            taxonomic_abundance: self.rel_abund.ok_or_else(|| {
                PyRuntimeError::new_err("profile result without taxonomic abundance")
            })?,
            kmers_reassigned: self.kmers_lost.ok_or_else(|| {
                PyRuntimeError::new_err("profile result without reassigned k-mers")
            })?,
        };
        let genome = Py::new(
            py,
            PyClassInitializer::from(Sketch).add_subclass(GenomeSketch::from(self.genome)),
        )?;
        Py::new(
            py,
            PyClassInitializer::from(AniResult {
                final_est_ani: self.final_est_ani,
                naive_ani: self.naive_ani,
                final_est_cov: self.final_est_cov,
                mean_cov: self.mean_cov,
                median_cov: self.median_cov,
                containment_index: self.containment_index,
                lambda: self.lambda,
                ani_ci: self.ani_ci,
                lambda_ci: self.lambda_ci,
                genome,
            })
            .add_subclass(profile),
        )
    }
}

/// Run the ``sylph query`` (or ``sylph profile``, if ``args.pseudotax``)
/// pipeline for a sample against a list of genome sketches.
///
/// Results are mapped back to the sketches they were computed from by
/// pointer identity, since several sketches can share the same name.
fn run_pipeline(
    args: &sylph::cmdline::ContainArgs,
    sketches: &[Arc<sylph::types::GenomeSketch>],
    sample: &sylph::types::SequencesSketch,
) -> Vec<Hit> {
    // estimate sample kmer identity (seq_id is a percentage, as in sylph)
    let kmer_id_opt = if let Some(x) = args.seq_id {
        Some((x / 100.0).powf(sample.k as f64))
    } else {
        sylph::contain::get_kmer_identity(sample, args.estimate_unknown)
    };

    // extract all matching kmers
    let mut stats = sketches
        .par_iter()
        .flat_map_iter(|sketch| sylph::contain::get_stats(args, sketch, sample, None, false))
        .collect::<Vec<_>>();

    // estimate true coverage
    sylph::contain::estimate_true_cov(
        &mut stats,
        kmer_id_opt,
        args.estimate_unknown,
        sample.mean_read_length,
        sample.k,
    );

    if args.pseudotax {
        // reassign k-mers
        let winner_map = sylph::contain::winner_table(&stats, false);
        let remaining_stats = stats
            .iter()
            .map(|x| x.genome_sketch)
            .flat_map(|genome_sketch| {
                sylph::contain::get_stats(args, genome_sketch, sample, Some(&winner_map), false)
            })
            .collect();
        stats = sylph::contain::derep_if_reassign_threshold(
            &stats,
            remaining_stats,
            args.redundant_ani,
            sample.k,
        );
        sylph::contain::estimate_true_cov(
            &mut stats,
            kmer_id_opt,
            args.estimate_unknown,
            sample.mean_read_length,
            sample.k,
        );

        // estimate unknown if needed
        let mut bases_explained = 1.;
        if args.estimate_unknown {
            bases_explained = sylph::contain::estimate_covered_bases(
                &stats,
                sample,
                sample.mean_read_length,
                sample.k,
            );
        }

        // compute coverage and relative abundances
        let total_cov = stats.iter().map(|x| x.final_est_cov).sum::<f64>();
        let total_seq_cov = stats
            .iter()
            .map(|x| x.final_est_cov * x.genome_sketch.gn_size as f64)
            .sum::<f64>();
        for result in stats.iter_mut() {
            result.rel_abund = Some(result.final_est_cov / total_cov * 100.0);
            result.seq_abund = Some(
                result.final_est_cov * result.genome_sketch.gn_size as f64 / total_seq_cov * 100.0
                    * bases_explained,
            );
        }

        // sort by relative abundance
        stats.sort_by(|x, y| {
            y.rel_abund
                .unwrap()
                .partial_cmp(&x.rel_abund.unwrap())
                .unwrap()
        });
    } else {
        // sort by ANI
        stats.sort_by(|x, y| y.final_est_ani.partial_cmp(&x.final_est_ani).unwrap());
    }

    // take ownership of the results
    let genomes = sketches
        .iter()
        .map(|arc| (Arc::as_ptr(arc), arc))
        .collect::<std::collections::HashMap<_, _>>();
    stats.iter().map(|r| Hit::new(r, &genomes)).collect()
}

/// A ``sylph`` profiler.
///
/// Arguments:
///     database (`~pysylph.Database` or `~pysylph.TwoStageDatabase`): The
///         database of genome sketches to profile samples against.
///     minimum_ani (`float`, *optional*): The minimum adjusted ANI to
///         consider (0-100). Defaults to 90 for `Profiler.query` and 95
///         for `Profiler.profile`.
///     seq_id (`float`, *optional*): The sequence identity (%) of the
///         reads, only used with ``estimate_unknown``. Overrides the
///         automatic detection.
///     estimate_unknown (`bool`): Estimate the true coverage and scale the
///         sequence abundance by the estimated unknown sequence percentage.
///     min_number_kmers (`int`): Discard genomes with fewer than this many
///         sampled k-mers. The default was lowered from 50 to 10 in
///         ``sylph`` v1.0.0.
///     min_contain (`int`): Discard genomes with fewer than this many
///         contained k-mers. Introduced in ``sylph`` v1.0.0.
///     screen_ani (`float`): The minimum adjusted ANI (0-100) for a genome
///         to pass the first-stage screen of a `TwoStageDatabase`. This is
///         deliberately permissive, and ignored for a `Database`.
///
/// Raises:
///     `TypeError`: When ``database`` is neither a `Database` nor a
///         `TwoStageDatabase`.
///     `ValueError`: When ``minimum_ani`` or ``screen_ani`` is out of range.
///         The `query` and `profile` methods also raise it when the sample
///         sketch is incompatible with the database (different ``k``,
///         or sample ``c`` greater than the genome ``c``).
#[pyclass(module = "pysylph.lib", frozen)]
pub struct Profiler {
    minimum_ani: Option<f64>,
    seq_id: Option<f64>,
    estimate_unknown: bool,
    min_number_kmers: usize,
    min_contain: usize,
    screen_ani: f64,
    database: ProfilerDatabase,
}

/// The database a profiler runs against.
enum ProfilerDatabase {
    Plain(Py<Database>),
    TwoStage(Py<TwoStageDatabase>),
}

#[pymethods]
impl Profiler {
    #[new]
    #[pyo3(signature = (database, *, minimum_ani = None, seq_id = None, estimate_unknown = false, min_number_kmers = 10, min_contain = 7, screen_ani = sylph::constants::SCREEN_MIN_ANI_DEFAULT))]
    pub fn __new__<'py>(
        database: &Bound<'py, PyAny>,
        minimum_ani: Option<f64>,
        seq_id: Option<f64>,
        estimate_unknown: bool,
        min_number_kmers: usize,
        min_contain: usize,
        screen_ani: f64,
    ) -> PyResult<Self> {
        let database = if let Ok(db) = database.cast::<Database>() {
            ProfilerDatabase::Plain(db.clone().unbind())
        } else if let Ok(db) = database.cast::<TwoStageDatabase>() {
            ProfilerDatabase::TwoStage(db.clone().unbind())
        } else {
            return Err(PyTypeError::new_err(
                "expected Database or TwoStageDatabase for database",
            ));
        };
        if screen_ani < 0.0 || screen_ani > 100.0 {
            return Err(PyValueError::new_err(format!(
                "invalid value for screen_ani: {}",
                screen_ani
            )));
        }
        if let Some(m) = minimum_ani {
            if m < 0.0 || m > 100.0 {
                return Err(PyValueError::new_err(format!(
                    "invalid value for minimum_ani: {}",
                    m
                )));
            }
        }
        Ok(Self {
            minimum_ani,
            database,
            estimate_unknown,
            seq_id,
            min_number_kmers,
            min_contain,
            screen_ani,
        })
    }

    /// `~pysylph.Database` or `~pysylph.TwoStageDatabase`: The database
    /// samples are profiled against.
    #[getter]
    fn database<'py>(&self, py: Python<'py>) -> Py<PyAny> {
        match &self.database {
            ProfilerDatabase::Plain(db) => db.clone_ref(py).into_any(),
            ProfilerDatabase::TwoStage(db) => db.clone_ref(py).into_any(),
        }
    }

    /// Run an ANI containment query for the given sample.
    ///
    /// Arguments:
    ///     sample (`~pysylph.SampleSketch`): The sketched sample to profile.
    ///
    /// Returns:
    ///     `list` of `~pysylph.AniResult`: The list of hits found in the
    ///     database for the sample.
    ///
    #[pyo3(signature = (sample))]
    fn query<'py>(&self, sample: PyRef<'py, SampleSketch>) -> PyResult<Vec<Py<AniResult>>> {
        let py = sample.py();
        let hits = self.run(py, false, &sample)?;
        hits.into_iter().map(|hit| hit.into_ani_result(py)).collect()
    }

    /// Build a taxonomic profile for the given sample.
    ///
    /// Arguments:
    ///     sample (`~pysylph.SampleSketch`): The sketched sample to profile.
    ///
    /// Returns:
    ///     `list` of `~pysylph.ProfileResult`: The list of hits found in
    ///     the database for the sample, with abundance estimates in
    ///     addition to coveraged-corrected ANI.
    ///
    #[pyo3(signature = (sample))]
    fn profile<'py>(&self, sample: PyRef<'py, SampleSketch>) -> PyResult<Vec<Py<ProfileResult>>> {
        let py = sample.py();
        let hits = self.run(py, true, &sample)?;
        hits.into_iter()
            .map(|hit| hit.into_profile_result(py))
            .collect()
    }
}

impl Profiler {
    fn run<'py>(
        &self,
        py: Python<'py>,
        pseudotax: bool,
        sample: &SampleSketch,
    ) -> PyResult<Vec<Hit>> {
        check_compatible(&self.database, sample)?;
        let mut args = sylph::cmdline::ContainArgs::from_cli_defaults(pseudotax);
        args.minimum_ani = self.minimum_ani;
        args.min_number_kmers = self.min_number_kmers as f64;
        args.min_contain = self.min_contain;
        args.estimate_unknown = self.estimate_unknown;
        args.seq_id = self.seq_id;
        args.no_sketch_pipeline = true;
        args.screen_ani = self.screen_ani;
        args.screen_dump = None;

        match &self.database {
            ProfilerDatabase::Plain(database) => {
                let sketches = &database.get().sketches;
                Ok(py.detach(|| run_pipeline(&args, sketches, &sample.sketch)))
            }
            ProfilerDatabase::TwoStage(database) => {
                let database = &database.get().db;
                Ok(py.detach(|| {
                    // screen the database, and run the pipeline on the survivors
                    let mut dense =
                        sylph::contain::compute_dense_survivors(&args, database, &sample.sketch);
                    // survivors are collected by parallel workers: sort them
                    // to get reproducible results
                    dense.sort_by(|a, b| {
                        (&a.file_name, &a.first_contig_name)
                            .cmp(&(&b.file_name, &b.first_contig_name))
                    });
                    let sketches = dense.into_iter().map(Arc::new).collect::<Vec<_>>();
                    run_pipeline(&args, &sketches, &sample.sketch)
                }))
            }
        }
    }
}

// --- Initializer -------------------------------------------------------------

/// PyO3 bindings to ``sylph``, an ultrafast taxonomic profiler.
#[pymodule]
#[pyo3(name = "lib")]
pub fn init(_py: Python, m: Bound<PyModule>) -> PyResult<()> {
    m.add("__package__", "pysylph")?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__author__", env!("CARGO_PKG_AUTHORS").replace(':', "\n"))?;

    m.add_class::<Sketch>()?;
    m.add_class::<GenomeSketch>()?;
    m.add_class::<SampleSketch>()?;

    m.add_class::<Database>()?;
    m.add_class::<DatabaseFile>()?;
    m.add_class::<TwoStageDatabase>()?;

    m.add_class::<Sketcher>()?;
    m.add_class::<Profiler>()?;

    m.add_class::<AniResult>()?;
    m.add_class::<ProfileResult>()?;

    Ok(())
}
