use lattice_core::{
    read_basis_from_file as core_read_basis_from_file,
    read_basis_from_str as core_read_basis_from_str,
    read_graph_from_file as core_read_graph_from_file,
    read_graph_from_str as core_read_graph_from_str,
    read_unitcell_from_file as core_read_unitcell_from_file,
    read_unitcell_from_str as core_read_unitcell_from_str,
    write_basis_to_string as core_write_basis_to_string,
    write_graph_to_string as core_write_graph_to_string,
    write_unitcell_to_string as core_write_unitcell_to_string, Basis, BasisMatrix, Boundary,
    CoordinateVector, ExtentVector, Graph, OffsetVector, Unitcell, XmlError,
};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::{PyIndexError, PyOSError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::PathBuf;

#[pyclass(name = "Boundary", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
enum PyBoundary {
    Open,
    Periodic,
}

impl From<PyBoundary> for Boundary {
    fn from(value: PyBoundary) -> Self {
        match value {
            PyBoundary::Open => Boundary::Open,
            PyBoundary::Periodic => Boundary::Periodic,
        }
    }
}

#[pyclass(name = "Basis")]
struct PyBasis {
    inner: Basis,
}

#[pymethods]
impl PyBasis {
    #[new]
    fn new(matrix: Vec<Vec<f64>>) -> PyResult<Self> {
        Ok(Self {
            inner: core_value(|| Ok(Basis::new(basis_matrix_from_rows(matrix)?)))?,
        })
    }

    #[staticmethod]
    fn simple(dim: usize) -> PyResult<Self> {
        Ok(Self {
            inner: core_value(|| Ok(Basis::simple(dim)))?,
        })
    }

    #[staticmethod]
    fn from_xml(xml: &str, name: &str) -> PyResult<Self> {
        read_basis_from_string(xml, name)
    }

    #[staticmethod]
    fn from_xml_file(path: PathBuf, name: &str) -> PyResult<Self> {
        read_basis_from_file(path, name)
    }

    fn to_xml(&self, name: &str) -> PyResult<String> {
        core_value(|| Ok(core_write_basis_to_string(name, &self.inner)))
    }

    #[getter]
    fn dimension(&self) -> usize {
        self.inner.dimension()
    }

    fn volume(&self) -> f64 {
        self.inner.volume()
    }

    fn matrix<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        PyArray2::from_vec2(py, &matrix_to_rows(self.inner.basis_vectors()))
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    fn matrix_list(&self) -> Vec<Vec<f64>> {
        matrix_to_rows(self.inner.basis_vectors())
    }

    fn __repr__(&self) -> String {
        format!("Basis(dimension={})", self.inner.dimension())
    }
}

#[pyclass(name = "Unitcell")]
struct PyUnitcell {
    inner: Unitcell,
}

#[pymethods]
impl PyUnitcell {
    #[new]
    fn new(dim: usize) -> PyResult<Self> {
        Ok(Self {
            inner: core_value(|| Ok(Unitcell::new(dim)))?,
        })
    }

    #[staticmethod]
    fn simple(dim: usize) -> PyResult<Self> {
        Ok(Self {
            inner: core_value(|| Ok(Unitcell::simple(dim)))?,
        })
    }

    #[staticmethod]
    fn from_xml(xml: &str, name: &str) -> PyResult<Self> {
        read_unitcell_from_string(xml, name)
    }

    #[staticmethod]
    fn from_xml_file(path: PathBuf, name: &str) -> PyResult<Self> {
        read_unitcell_from_file(path, name)
    }

    fn to_xml(&self, name: &str) -> PyResult<String> {
        core_value(|| Ok(core_write_unitcell_to_string(name, &self.inner)))
    }

    #[getter]
    fn dimension(&self) -> usize {
        self.inner.dimension()
    }

    #[getter]
    fn num_sites(&self) -> usize {
        self.inner.num_sites()
    }

    #[getter]
    fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }

    fn max_neighbors(&self) -> usize {
        self.inner.max_neighbors()
    }

    #[pyo3(signature = (coordinate, site_type=0))]
    fn add_site(&mut self, coordinate: Vec<f64>, site_type: i32) -> PyResult<usize> {
        validate_vector_len(coordinate.len(), self.inner.dimension(), "site coordinate")?;
        if coordinate.iter().any(|value| *value < 0.0 || *value >= 1.0) {
            return Err(PyValueError::new_err(
                "site coordinate values must be in [0, 1)",
            ));
        }
        core_value(|| {
            Ok(self
                .inner
                .add_site(CoordinateVector::from_vec(coordinate), site_type))
        })
    }

    #[pyo3(signature = (source, target, target_offset, bond_type=0))]
    fn add_bond(
        &mut self,
        source: usize,
        target: usize,
        target_offset: Vec<i64>,
        bond_type: i32,
    ) -> PyResult<usize> {
        validate_site_index(source, self.inner.num_sites(), "source")?;
        validate_site_index(target, self.inner.num_sites(), "target")?;
        validate_vector_len(target_offset.len(), self.inner.dimension(), "target offset")?;
        core_value(|| {
            Ok(self.inner.add_bond(
                source,
                target,
                OffsetVector::from_vec(target_offset),
                bond_type,
            ))
        })
    }

    fn site_type(&self, site: usize) -> PyResult<i32> {
        validate_site_index(site, self.inner.num_sites(), "site")?;
        Ok(self.inner.site(site).site_type)
    }

    fn coordinate(&self, site: usize) -> PyResult<Vec<f64>> {
        validate_site_index(site, self.inner.num_sites(), "site")?;
        Ok(self.inner.site(site).coordinate.iter().copied().collect())
    }

    fn bond_source(&self, bond: usize) -> PyResult<usize> {
        validate_bond_index(bond, self.inner.num_bonds())?;
        Ok(self.inner.bond(bond).source)
    }

    fn bond_target(&self, bond: usize) -> PyResult<usize> {
        validate_bond_index(bond, self.inner.num_bonds())?;
        Ok(self.inner.bond(bond).target)
    }

    fn bond_type(&self, bond: usize) -> PyResult<i32> {
        validate_bond_index(bond, self.inner.num_bonds())?;
        Ok(self.inner.bond(bond).bond_type)
    }

    fn bond_offset(&self, bond: usize) -> PyResult<Vec<i64>> {
        validate_bond_index(bond, self.inner.num_bonds())?;
        Ok(self
            .inner
            .bond(bond)
            .target_offset
            .iter()
            .copied()
            .collect())
    }

    fn __repr__(&self) -> String {
        format!(
            "Unitcell(dimension={}, num_sites={}, num_bonds={})",
            self.inner.dimension(),
            self.inner.num_sites(),
            self.inner.num_bonds()
        )
    }
}

#[pyclass(name = "Graph")]
struct PyGraph {
    inner: Graph,
}

#[pymethods]
impl PyGraph {
    #[new]
    fn new(dim: usize) -> PyResult<Self> {
        Ok(Self {
            inner: core_value(|| Ok(Graph::new(dim)))?,
        })
    }

    #[staticmethod]
    fn simple(dim: usize, length: usize) -> PyResult<Self> {
        if length == 0 {
            return Err(PyValueError::new_err("length must be positive"));
        }
        Ok(Self {
            inner: core_value(|| Ok(Graph::simple(dim, length)))?,
        })
    }

    #[staticmethod]
    fn fully_connected(num_sites: usize) -> PyResult<Self> {
        Ok(Self {
            inner: core_value(|| Ok(Graph::fully_connected(num_sites)))?,
        })
    }

    #[staticmethod]
    fn from_xml(xml: &str, name: &str) -> PyResult<Self> {
        read_graph_from_string(xml, name)
    }

    #[staticmethod]
    fn from_xml_file(path: PathBuf, name: &str) -> PyResult<Self> {
        read_graph_from_file(path, name)
    }

    fn to_xml(&self, name: &str) -> PyResult<String> {
        core_value(|| Ok(core_write_graph_to_string(name, &self.inner)))
    }

    #[staticmethod]
    fn from_basis_unitcell_extent(
        basis: PyRef<'_, PyBasis>,
        unitcell: PyRef<'_, PyUnitcell>,
        extent: Vec<i64>,
        boundary: Vec<PyBoundary>,
    ) -> PyResult<Self> {
        validate_lattice_inputs(&basis.inner, &unitcell.inner, &extent, boundary.len())?;
        let extent = ExtentVector::from_vec(extent);
        let boundary: Vec<Boundary> = boundary.into_iter().map(Boundary::from).collect();
        Ok(Self {
            inner: core_value(|| {
                Ok(Graph::from_basis_unitcell_extent(
                    &basis.inner,
                    &unitcell.inner,
                    &extent,
                    &boundary,
                ))
            })?,
        })
    }

    #[staticmethod]
    fn from_basis_unitcell_length(
        basis: PyRef<'_, PyBasis>,
        unitcell: PyRef<'_, PyUnitcell>,
        length: usize,
        boundary: PyBoundary,
    ) -> PyResult<Self> {
        if basis.inner.dimension() != unitcell.inner.dimension() {
            return Err(PyValueError::new_err(
                "basis and unitcell dimensions must match",
            ));
        }
        if length == 0 {
            return Err(PyValueError::new_err("length must be positive"));
        }
        Ok(Self {
            inner: core_value(|| {
                Ok(Graph::from_basis_unitcell_length(
                    &basis.inner,
                    &unitcell.inner,
                    length,
                    Boundary::from(boundary),
                ))
            })?,
        })
    }

    #[getter]
    fn dimension(&self) -> usize {
        self.inner.dimension()
    }

    #[getter]
    fn num_sites(&self) -> usize {
        self.inner.num_sites()
    }

    #[getter]
    fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }

    fn site_type(&self, site: usize) -> PyResult<i32> {
        validate_site_index(site, self.inner.num_sites(), "site")?;
        Ok(self.inner.site_type(site))
    }

    fn coordinate(&self, site: usize) -> PyResult<Vec<f64>> {
        validate_site_index(site, self.inner.num_sites(), "site")?;
        Ok(self.inner.coordinate(site).iter().copied().collect())
    }

    fn num_neighbors(&self, site: usize) -> PyResult<usize> {
        validate_site_index(site, self.inner.num_sites(), "site")?;
        Ok(self.inner.num_neighbors(site))
    }

    fn neighbor(&self, site: usize, neighbor: usize) -> PyResult<usize> {
        validate_site_index(site, self.inner.num_sites(), "site")?;
        validate_neighbor_index(neighbor, self.inner.num_neighbors(site))?;
        Ok(self.inner.neighbor(site, neighbor))
    }

    fn neighbor_bond(&self, site: usize, neighbor: usize) -> PyResult<usize> {
        validate_site_index(site, self.inner.num_sites(), "site")?;
        validate_neighbor_index(neighbor, self.inner.num_neighbors(site))?;
        Ok(self.inner.neighbor_bond(site, neighbor))
    }

    fn bond_type(&self, bond: usize) -> PyResult<i32> {
        validate_bond_index(bond, self.inner.num_bonds())?;
        Ok(self.inner.bond_type(bond))
    }

    fn source(&self, bond: usize) -> PyResult<usize> {
        validate_bond_index(bond, self.inner.num_bonds())?;
        Ok(self.inner.source(bond))
    }

    fn target(&self, bond: usize) -> PyResult<usize> {
        validate_bond_index(bond, self.inner.num_bonds())?;
        Ok(self.inner.target(bond))
    }

    fn edge_sites(&self, bond: usize) -> PyResult<(usize, usize)> {
        validate_bond_index(bond, self.inner.num_bonds())?;
        Ok(self.inner.edge_sites(bond))
    }

    #[pyo3(signature = (coordinate, site_type=0))]
    fn add_site(&mut self, coordinate: Vec<f64>, site_type: i32) -> PyResult<usize> {
        validate_vector_len(coordinate.len(), self.inner.dimension(), "site coordinate")?;
        core_value(|| {
            Ok(self
                .inner
                .add_site(CoordinateVector::from_vec(coordinate), site_type))
        })
    }

    #[pyo3(signature = (source, target, bond_type=0))]
    fn add_bond(&mut self, source: usize, target: usize, bond_type: i32) -> PyResult<usize> {
        validate_site_index(source, self.inner.num_sites(), "source")?;
        validate_site_index(target, self.inner.num_sites(), "target")?;
        if source == target {
            return Err(PyValueError::new_err("self loop is not allowed"));
        }
        core_value(|| Ok(self.inner.add_bond(source, target, bond_type)))
    }

    fn coordinates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        PyArray2::from_vec2(py, &self.coordinates_vec())
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    fn site_types<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i32>> {
        self.site_types_vec().into_pyarray(py)
    }

    fn edges<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<usize>>> {
        PyArray2::from_vec2(py, &self.edges_vec())
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    fn bond_types<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i32>> {
        self.bond_types_vec().into_pyarray(py)
    }

    fn coordinates_list(&self) -> Vec<Vec<f64>> {
        self.coordinates_vec()
    }

    fn site_types_list(&self) -> Vec<i32> {
        self.site_types_vec()
    }

    fn edges_list(&self) -> Vec<Vec<usize>> {
        self.edges_vec()
    }

    fn bond_types_list(&self) -> Vec<i32> {
        self.bond_types_vec()
    }

    fn __repr__(&self) -> String {
        format!(
            "Graph(dimension={}, num_sites={}, num_bonds={})",
            self.inner.dimension(),
            self.inner.num_sites(),
            self.inner.num_bonds()
        )
    }
}

impl PyGraph {
    fn coordinates_vec(&self) -> Vec<Vec<f64>> {
        (0..self.inner.num_sites())
            .map(|site| self.inner.coordinate(site).iter().copied().collect())
            .collect()
    }

    fn site_types_vec(&self) -> Vec<i32> {
        (0..self.inner.num_sites())
            .map(|site| self.inner.site_type(site))
            .collect()
    }

    fn edges_vec(&self) -> Vec<Vec<usize>> {
        (0..self.inner.num_bonds())
            .map(|bond| {
                let (source, target) = self.inner.edge_sites(bond);
                vec![source, target]
            })
            .collect()
    }

    fn bond_types_vec(&self) -> Vec<i32> {
        (0..self.inner.num_bonds())
            .map(|bond| self.inner.bond_type(bond))
            .collect()
    }
}

#[pyfunction]
fn read_basis_from_string(xml: &str, name: &str) -> PyResult<PyBasis> {
    Ok(PyBasis {
        inner: core_result(|| core_read_basis_from_str(xml, name))?,
    })
}

#[pyfunction]
fn read_unitcell_from_string(xml: &str, name: &str) -> PyResult<PyUnitcell> {
    Ok(PyUnitcell {
        inner: core_result(|| core_read_unitcell_from_str(xml, name))?,
    })
}

#[pyfunction]
fn read_graph_from_string(xml: &str, name: &str) -> PyResult<PyGraph> {
    Ok(PyGraph {
        inner: core_result(|| core_read_graph_from_str(xml, name))?,
    })
}

#[pyfunction]
fn read_basis_from_file(path: PathBuf, name: &str) -> PyResult<PyBasis> {
    Ok(PyBasis {
        inner: core_result(|| core_read_basis_from_file(path, name))?,
    })
}

#[pyfunction]
fn read_unitcell_from_file(path: PathBuf, name: &str) -> PyResult<PyUnitcell> {
    Ok(PyUnitcell {
        inner: core_result(|| core_read_unitcell_from_file(path, name))?,
    })
}

#[pyfunction]
fn read_graph_from_file(path: PathBuf, name: &str) -> PyResult<PyGraph> {
    Ok(PyGraph {
        inner: core_result(|| core_read_graph_from_file(path, name))?,
    })
}

#[pyfunction]
fn write_basis_to_string(name: &str, basis: PyRef<'_, PyBasis>) -> PyResult<String> {
    core_value(|| Ok(core_write_basis_to_string(name, &basis.inner)))
}

#[pyfunction]
fn write_unitcell_to_string(name: &str, unitcell: PyRef<'_, PyUnitcell>) -> PyResult<String> {
    core_value(|| Ok(core_write_unitcell_to_string(name, &unitcell.inner)))
}

#[pyfunction]
fn write_graph_to_string(name: &str, graph: PyRef<'_, PyGraph>) -> PyResult<String> {
    core_value(|| Ok(core_write_graph_to_string(name, &graph.inner)))
}

#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[pymodule]
fn lattice(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<PyBoundary>()?;
    m.add_class::<PyBasis>()?;
    m.add_class::<PyUnitcell>()?;
    m.add_class::<PyGraph>()?;
    m.add_function(wrap_pyfunction!(read_basis_from_string, m)?)?;
    m.add_function(wrap_pyfunction!(read_unitcell_from_string, m)?)?;
    m.add_function(wrap_pyfunction!(read_graph_from_string, m)?)?;
    m.add_function(wrap_pyfunction!(read_basis_from_file, m)?)?;
    m.add_function(wrap_pyfunction!(read_unitcell_from_file, m)?)?;
    m.add_function(wrap_pyfunction!(read_graph_from_file, m)?)?;
    m.add_function(wrap_pyfunction!(write_basis_to_string, m)?)?;
    m.add_function(wrap_pyfunction!(write_unitcell_to_string, m)?)?;
    m.add_function(wrap_pyfunction!(write_graph_to_string, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    Ok(())
}

fn basis_matrix_from_rows(rows: Vec<Vec<f64>>) -> PyResult<BasisMatrix> {
    let dim = rows.len();
    if rows.iter().any(|row| row.len() != dim) {
        return Err(PyValueError::new_err("basis matrix must be square"));
    }
    let values: Vec<f64> = rows.into_iter().flatten().collect();
    Ok(BasisMatrix::from_row_slice(dim, dim, &values))
}

fn matrix_to_rows(matrix: &BasisMatrix) -> Vec<Vec<f64>> {
    (0..matrix.nrows())
        .map(|row| {
            (0..matrix.ncols())
                .map(|column| matrix[(row, column)])
                .collect()
        })
        .collect()
}

fn validate_lattice_inputs(
    basis: &Basis,
    unitcell: &Unitcell,
    extent: &[i64],
    boundary_len: usize,
) -> PyResult<()> {
    let dim = unitcell.dimension();
    if basis.dimension() != dim {
        return Err(PyValueError::new_err(
            "basis and unitcell dimensions must match",
        ));
    }
    validate_vector_len(extent.len(), dim, "extent")?;
    if boundary_len != dim {
        return Err(PyValueError::new_err(
            "boundary dimension must match unitcell dimension",
        ));
    }
    if extent.iter().any(|value| *value <= 0) {
        return Err(PyValueError::new_err("extent values must be positive"));
    }
    Ok(())
}

fn validate_vector_len(actual: usize, expected: usize, name: &str) -> PyResult<()> {
    if actual != expected {
        return Err(PyValueError::new_err(format!(
            "{name} dimension mismatch: expected {expected}, got {actual}"
        )));
    }
    Ok(())
}

fn validate_site_index(index: usize, len: usize, name: &str) -> PyResult<()> {
    if index >= len {
        return Err(PyIndexError::new_err(format!("{name} index out of range")));
    }
    Ok(())
}

fn validate_bond_index(index: usize, len: usize) -> PyResult<()> {
    if index >= len {
        return Err(PyIndexError::new_err("bond index out of range"));
    }
    Ok(())
}

fn validate_neighbor_index(index: usize, len: usize) -> PyResult<()> {
    if index >= len {
        return Err(PyIndexError::new_err("neighbor index out of range"));
    }
    Ok(())
}

fn xml_error_to_py_err(error: XmlError) -> PyErr {
    let message = error.to_string();
    match error {
        XmlError::Io(_) => PyOSError::new_err(message),
        _ => PyValueError::new_err(message),
    }
}

fn core_result<T>(f: impl FnOnce() -> Result<T, XmlError>) -> PyResult<T> {
    core_value(|| f().map_err(xml_error_to_py_err))
}

fn core_value<T>(f: impl FnOnce() -> PyResult<T>) -> PyResult<T> {
    match catch_unwind(AssertUnwindSafe(f)) {
        Ok(result) => result,
        Err(_) => Err(PyRuntimeError::new_err("lattice-core panicked")),
    }
}
