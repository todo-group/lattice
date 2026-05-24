use lattice_core::{read_basis_from_str, read_graph_from_str, read_unitcell_from_str, write_basis_to_string, write_graph_to_string, write_unitcell_to_string, Basis, BasisMatrix, CoordinateVector, OffsetVector, Graph, Unitcell};
use std::cell::RefCell;
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_double, c_int};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::ptr;

thread_local! {
  static LAST_ERROR: RefCell<Option<String>> = const { RefCell::new(None) };
}

fn set_last_error(message: impl Into<String>) {
  LAST_ERROR.with(|slot| {
    *slot.borrow_mut() = Some(message.into());
  });
}

fn clear_last_error() {
  LAST_ERROR.with(|slot| {
    *slot.borrow_mut() = None;
  });
}

fn copy_last_error() -> Option<String> {
  LAST_ERROR.with(|slot| slot.borrow().clone())
}

#[repr(C)]
pub struct lattice_basis_raw {
  pub dim: usize,
  pub values_len: usize,
  pub values: *mut c_double,
}

#[repr(C)]
pub struct lattice_unitcell_raw {
  pub dim: usize,
  pub num_sites: usize,
  pub site_types: *mut c_int,
  pub site_coordinates_len: usize,
  pub site_coordinates: *mut c_double,
  pub num_bonds: usize,
  pub bond_sources: *mut usize,
  pub bond_targets: *mut usize,
  pub bond_types: *mut c_int,
  pub bond_offsets_len: usize,
  pub bond_offsets: *mut i64,
}

#[repr(C)]
pub struct lattice_graph_raw {
  pub dim: usize,
  pub num_sites: usize,
  pub site_types: *mut c_int,
  pub site_coordinates_len: usize,
  pub site_coordinates: *mut c_double,
  pub num_bonds: usize,
  pub bond_sources: *mut usize,
  pub bond_targets: *mut usize,
  pub bond_types: *mut c_int,
}

#[no_mangle]
pub extern "C" fn lattice_string_free(ptr: *mut c_char) {
  ffi_void(|| {
    if !ptr.is_null() {
      unsafe {
        drop(CString::from_raw(ptr));
      }
    }
  });
}

#[no_mangle]
pub extern "C" fn lattice_last_error_message() -> *mut c_char {
  ffi_ptr("lattice_last_error_message", || {
    let Some(message) = copy_last_error() else {
      return ptr::null_mut();
    };
    string_to_c(message)
  })
}

#[no_mangle]
pub extern "C" fn lattice_basis_raw_free(ptr: *mut lattice_basis_raw) {
  ffi_void(|| {
    if ptr.is_null() {
      return;
    }
    unsafe {
      let raw = Box::from_raw(ptr);
      if !raw.values.is_null() {
        drop(Vec::from_raw_parts(raw.values, raw.values_len, raw.values_len));
      }
    }
  });
}

#[no_mangle]
pub extern "C" fn lattice_unitcell_raw_free(ptr: *mut lattice_unitcell_raw) {
  ffi_void(|| {
    if ptr.is_null() {
      return;
    }
    unsafe {
      let raw = Box::from_raw(ptr);
      if !raw.site_types.is_null() {
        drop(Vec::from_raw_parts(raw.site_types, raw.num_sites, raw.num_sites));
      }
      if !raw.site_coordinates.is_null() {
        drop(Vec::from_raw_parts(raw.site_coordinates, raw.site_coordinates_len, raw.site_coordinates_len));
      }
      if !raw.bond_sources.is_null() {
        drop(Vec::from_raw_parts(raw.bond_sources, raw.num_bonds, raw.num_bonds));
      }
      if !raw.bond_targets.is_null() {
        drop(Vec::from_raw_parts(raw.bond_targets, raw.num_bonds, raw.num_bonds));
      }
      if !raw.bond_types.is_null() {
        drop(Vec::from_raw_parts(raw.bond_types, raw.num_bonds, raw.num_bonds));
      }
      if !raw.bond_offsets.is_null() {
        drop(Vec::from_raw_parts(raw.bond_offsets, raw.bond_offsets_len, raw.bond_offsets_len));
      }
    }
  });
}

#[no_mangle]
pub extern "C" fn lattice_graph_raw_free(ptr: *mut lattice_graph_raw) {
  ffi_void(|| {
    if ptr.is_null() {
      return;
    }
    unsafe {
      let raw = Box::from_raw(ptr);
      if !raw.site_types.is_null() {
        drop(Vec::from_raw_parts(raw.site_types, raw.num_sites, raw.num_sites));
      }
      if !raw.site_coordinates.is_null() {
        drop(Vec::from_raw_parts(raw.site_coordinates, raw.site_coordinates_len, raw.site_coordinates_len));
      }
      if !raw.bond_sources.is_null() {
        drop(Vec::from_raw_parts(raw.bond_sources, raw.num_bonds, raw.num_bonds));
      }
      if !raw.bond_targets.is_null() {
        drop(Vec::from_raw_parts(raw.bond_targets, raw.num_bonds, raw.num_bonds));
      }
      if !raw.bond_types.is_null() {
        drop(Vec::from_raw_parts(raw.bond_types, raw.num_bonds, raw.num_bonds));
      }
    }
  });
}

#[no_mangle]
pub extern "C" fn lattice_basis_from_xml(xml: *const c_char, name: *const c_char) -> *mut lattice_basis_raw {
  ffi_ptr("lattice_basis_from_xml", || {
    with_strs(xml, name, |xml, name| {
      let basis = read_basis_from_str(xml, name).ok()?;
      Some(box_basis_raw(&basis))
    })
  })
}

#[no_mangle]
pub extern "C" fn lattice_unitcell_from_xml(xml: *const c_char, name: *const c_char) -> *mut lattice_unitcell_raw {
  ffi_ptr("lattice_unitcell_from_xml", || {
    with_strs(xml, name, |xml, name| {
      let cell = read_unitcell_from_str(xml, name).ok()?;
      Some(box_unitcell_raw(&cell))
    })
  })
}

#[no_mangle]
pub extern "C" fn lattice_graph_from_xml(xml: *const c_char, name: *const c_char) -> *mut lattice_graph_raw {
  ffi_ptr("lattice_graph_from_xml", || {
    with_strs(xml, name, |xml, name| {
      let graph = read_graph_from_str(xml, name).ok()?;
      Some(box_graph_raw(&graph))
    })
  })
}

#[no_mangle]
pub extern "C" fn lattice_basis_to_xml(name: *const c_char, raw: *const lattice_basis_raw) -> *mut c_char {
  ffi_ptr("lattice_basis_to_xml", || unsafe {
    if name.is_null() {
      return ptr::null_mut();
    }
    let name = match CStr::from_ptr(name).to_str() {
      Ok(value) => value,
      Err(_) => return ptr::null_mut(),
    };
    let raw = match raw.as_ref() {
      Some(value) => value,
      None => return ptr::null_mut(),
    };
    let basis = match basis_from_raw(raw) {
      Some(value) => value,
      None => return ptr::null_mut(),
    };
    string_to_c(write_basis_to_string(name, &basis))
  })
}

#[no_mangle]
pub extern "C" fn lattice_unitcell_to_xml(name: *const c_char, raw: *const lattice_unitcell_raw) -> *mut c_char {
  ffi_ptr("lattice_unitcell_to_xml", || unsafe {
    if name.is_null() {
      return ptr::null_mut();
    }
    let name = match CStr::from_ptr(name).to_str() {
      Ok(value) => value,
      Err(_) => return ptr::null_mut(),
    };
    let raw = match raw.as_ref() {
      Some(value) => value,
      None => return ptr::null_mut(),
    };
    let cell = match unitcell_from_raw(raw) {
      Some(value) => value,
      None => return ptr::null_mut(),
    };
    string_to_c(write_unitcell_to_string(name, &cell))
  })
}

#[no_mangle]
pub extern "C" fn lattice_graph_to_xml(name: *const c_char, raw: *const lattice_graph_raw) -> *mut c_char {
  ffi_ptr("lattice_graph_to_xml", || unsafe {
    if name.is_null() {
      return ptr::null_mut();
    }
    let name = match CStr::from_ptr(name).to_str() {
      Ok(value) => value,
      Err(_) => return ptr::null_mut(),
    };
    let raw = match raw.as_ref() {
      Some(value) => value,
      None => return ptr::null_mut(),
    };
    let graph = match graph_from_raw(raw) {
      Some(value) => value,
      None => return ptr::null_mut(),
    };
    string_to_c(write_graph_to_string(name, &graph))
  })
}

fn ffi_ptr<T>(api_name: &str, f: impl FnOnce() -> *mut T) -> *mut T {
  match catch_unwind(AssertUnwindSafe(f)) {
    Ok(ptr) => {
      if ptr.is_null() {
        set_last_error(format!("{api_name} failed"));
      } else {
        clear_last_error();
      }
      ptr
    }
    Err(_) => {
      set_last_error(format!("{api_name} panicked"));
      ptr::null_mut()
    }
  }
}

fn ffi_void(f: impl FnOnce()) {
  let _ = catch_unwind(AssertUnwindSafe(f));
}

fn with_strs<T>(xml: *const c_char, name: *const c_char, f: impl FnOnce(&str, &str) -> Option<T>) -> *mut T {
  unsafe {
    if xml.is_null() || name.is_null() {
      return ptr::null_mut();
    }
    let xml = match CStr::from_ptr(xml).to_str() {
      Ok(value) => value,
      Err(_) => return ptr::null_mut(),
    };
    let name = match CStr::from_ptr(name).to_str() {
      Ok(value) => value,
      Err(_) => return ptr::null_mut(),
    };
    match f(xml, name) {
      Some(value) => Box::into_raw(Box::new(value)),
      None => ptr::null_mut(),
    }
  }
}

fn string_to_c(text: String) -> *mut c_char {
  CString::new(text).ok().map_or(ptr::null_mut(), CString::into_raw)
}

fn boxed_slice<T>(values: Vec<T>) -> (*mut T, usize) {
  let mut values = values.into_boxed_slice();
  let len = values.len();
  let ptr = values.as_mut_ptr();
  std::mem::forget(values);
  (ptr, len)
}

fn box_basis_raw(basis: &Basis) -> lattice_basis_raw {
  let mut values = Vec::with_capacity(basis.dimension() * basis.dimension());
  for row in 0..basis.dimension() {
    for column in 0..basis.dimension() {
      values.push(basis.basis_vectors()[(row, column)]);
    }
  }
  let (values, values_len) = boxed_slice(values);
  lattice_basis_raw { dim: basis.dimension(), values_len, values }
}

fn basis_from_raw(raw: &lattice_basis_raw) -> Option<Basis> {
  if raw.dim == 0 || raw.values_len != raw.dim * raw.dim || raw.values.is_null() {
    return None;
  }
  unsafe {
    let values = std::slice::from_raw_parts(raw.values, raw.values_len);
    let mut matrix = BasisMatrix::zeros(raw.dim, raw.dim);
    for row in 0..raw.dim {
      for column in 0..raw.dim {
        matrix[(row, column)] = values[row * raw.dim + column];
      }
    }
    Some(Basis::new(matrix))
  }
}

fn box_unitcell_raw(cell: &Unitcell) -> lattice_unitcell_raw {
  let mut site_types = Vec::with_capacity(cell.num_sites());
  let mut site_coordinates = Vec::with_capacity(cell.num_sites() * cell.dimension());
  for site in 0..cell.num_sites() {
    site_types.push(cell.site(site).site_type as c_int);
    for value in cell.site(site).coordinate.iter() {
      site_coordinates.push(*value);
    }
  }
  let mut bond_sources = Vec::with_capacity(cell.num_bonds());
  let mut bond_targets = Vec::with_capacity(cell.num_bonds());
  let mut bond_types = Vec::with_capacity(cell.num_bonds());
  let mut bond_offsets = Vec::with_capacity(cell.num_bonds() * cell.dimension());
  for bond in 0..cell.num_bonds() {
    let bond = cell.bond(bond);
    bond_sources.push(bond.source);
    bond_targets.push(bond.target);
    bond_types.push(bond.bond_type as c_int);
    for value in bond.target_offset.iter() {
      bond_offsets.push(*value);
    }
  }
  let (site_types, _) = boxed_slice(site_types);
  let (site_coordinates, site_coordinates_len) = boxed_slice(site_coordinates);
  let (bond_sources, _) = boxed_slice(bond_sources);
  let (bond_targets, _) = boxed_slice(bond_targets);
  let (bond_types, _) = boxed_slice(bond_types);
  let (bond_offsets, bond_offsets_len) = boxed_slice(bond_offsets);
  lattice_unitcell_raw {
    dim: cell.dimension(),
    num_sites: cell.num_sites(),
    site_types,
    site_coordinates_len,
    site_coordinates,
    num_bonds: cell.num_bonds(),
    bond_sources,
    bond_targets,
    bond_types,
    bond_offsets_len,
    bond_offsets,
  }
}

fn unitcell_from_raw(raw: &lattice_unitcell_raw) -> Option<Unitcell> {
  if raw.site_types.is_null() || raw.site_coordinates.is_null() || raw.bond_sources.is_null() || raw.bond_targets.is_null() || raw.bond_types.is_null() || raw.bond_offsets.is_null() {
    return None;
  }
  if raw.site_coordinates_len != raw.num_sites * raw.dim || raw.bond_offsets_len != raw.num_bonds * raw.dim {
    return None;
  }
  unsafe {
    let site_types = std::slice::from_raw_parts(raw.site_types, raw.num_sites);
    let site_coordinates = std::slice::from_raw_parts(raw.site_coordinates, raw.site_coordinates_len);
    let bond_sources = std::slice::from_raw_parts(raw.bond_sources, raw.num_bonds);
    let bond_targets = std::slice::from_raw_parts(raw.bond_targets, raw.num_bonds);
    let bond_types = std::slice::from_raw_parts(raw.bond_types, raw.num_bonds);
    let bond_offsets = std::slice::from_raw_parts(raw.bond_offsets, raw.bond_offsets_len);
    let mut cell = Unitcell::new(raw.dim);
    for site in 0..raw.num_sites {
      let mut coordinate = CoordinateVector::zeros(raw.dim);
      for axis in 0..raw.dim {
        coordinate[axis] = site_coordinates[site * raw.dim + axis];
      }
      cell.add_site(coordinate, site_types[site]);
    }
    for bond in 0..raw.num_bonds {
      let mut offset = OffsetVector::zeros(raw.dim);
      for axis in 0..raw.dim {
        offset[axis] = bond_offsets[bond * raw.dim + axis];
      }
      cell.add_bond(bond_sources[bond], bond_targets[bond], offset, bond_types[bond]);
    }
    Some(cell)
  }
}

fn box_graph_raw(graph: &Graph) -> lattice_graph_raw {
  let mut site_types = Vec::with_capacity(graph.num_sites());
  let mut site_coordinates = Vec::with_capacity(graph.num_sites() * graph.dimension());
  for site in 0..graph.num_sites() {
    site_types.push(graph.site_type(site) as c_int);
    for value in graph.coordinate(site).iter() {
      site_coordinates.push(*value);
    }
  }
  let mut bond_sources = Vec::with_capacity(graph.num_bonds());
  let mut bond_targets = Vec::with_capacity(graph.num_bonds());
  let mut bond_types = Vec::with_capacity(graph.num_bonds());
  for bond in 0..graph.num_bonds() {
    bond_sources.push(graph.source(bond));
    bond_targets.push(graph.target(bond));
    bond_types.push(graph.bond_type(bond) as c_int);
  }
  let (site_types, _) = boxed_slice(site_types);
  let (site_coordinates, site_coordinates_len) = boxed_slice(site_coordinates);
  let (bond_sources, _) = boxed_slice(bond_sources);
  let (bond_targets, _) = boxed_slice(bond_targets);
  let (bond_types, _) = boxed_slice(bond_types);
  lattice_graph_raw {
    dim: graph.dimension(),
    num_sites: graph.num_sites(),
    site_types,
    site_coordinates_len,
    site_coordinates,
    num_bonds: graph.num_bonds(),
    bond_sources,
    bond_targets,
    bond_types,
  }
}

fn graph_from_raw(raw: &lattice_graph_raw) -> Option<Graph> {
  if raw.site_types.is_null() || raw.site_coordinates.is_null() || raw.bond_sources.is_null() || raw.bond_targets.is_null() || raw.bond_types.is_null() {
    return None;
  }
  if raw.site_coordinates_len != raw.num_sites * raw.dim {
    return None;
  }
  unsafe {
    let site_types = std::slice::from_raw_parts(raw.site_types, raw.num_sites);
    let site_coordinates = std::slice::from_raw_parts(raw.site_coordinates, raw.site_coordinates_len);
    let bond_sources = std::slice::from_raw_parts(raw.bond_sources, raw.num_bonds);
    let bond_targets = std::slice::from_raw_parts(raw.bond_targets, raw.num_bonds);
    let bond_types = std::slice::from_raw_parts(raw.bond_types, raw.num_bonds);
    let mut graph = Graph::new(raw.dim);
    for site in 0..raw.num_sites {
      let mut coordinate = CoordinateVector::zeros(raw.dim);
      for axis in 0..raw.dim {
        coordinate[axis] = site_coordinates[site * raw.dim + axis];
      }
      graph.add_site(coordinate, site_types[site]);
    }
    for bond in 0..raw.num_bonds {
      graph.add_bond(bond_sources[bond], bond_targets[bond], bond_types[bond]);
    }
    Some(graph)
  }
}
