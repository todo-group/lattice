use crate::basis::Basis;
use crate::types::{Boundary, CoordinateVector, ExtentVector, OffsetVector};
use crate::unitcell::Unitcell;

#[derive(Clone, Debug, PartialEq)]
pub struct Site {
  pub site_type: i32,
  pub neighbors: Vec<usize>,
  pub neighbor_bonds: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Bond {
  pub source: usize,
  pub target: usize,
  pub bond_type: i32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Graph {
  dim: usize,
  sites: Vec<Site>,
  coordinates: Vec<CoordinateVector>,
  bonds: Vec<Bond>,
}

impl Graph {
  pub fn new(dim: usize) -> Self {
    Self { dim, sites: Vec::new(), coordinates: Vec::new(), bonds: Vec::new() }
  }

  pub fn simple(dim: usize, length: usize) -> Self {
    let basis = Basis::simple(dim);
    let cell = Unitcell::simple(dim);
    let extent = ExtentVector::from_element(dim, length as i64);
    let boundary = vec![Boundary::Periodic; dim];
    Self::from_basis_unitcell_extent(&basis, &cell, &extent, &boundary)
  }

  pub fn fully_connected(num_sites: usize) -> Self {
    let mut graph = Self::new(0);
    let pos = CoordinateVector::zeros(0);
    for _ in 0..num_sites {
      graph.add_site(pos.clone(), 0);
    }
    for source in 0..num_sites {
      for target in (source + 1)..num_sites {
        graph.add_bond(source, target, 0);
      }
    }
    graph
  }

  pub fn from_basis_unitcell_extent(
    basis: &Basis,
    cell: &Unitcell,
    extent: &ExtentVector,
    boundary: &[Boundary],
  ) -> Self {
    assert_eq!(cell.dimension(), basis.dimension(), "dimension mismatch");
    assert_eq!(cell.dimension(), extent.len(), "dimension mismatch");
    assert_eq!(cell.dimension(), boundary.len(), "dimension mismatch");
    for value in extent.iter() {
      assert!(*value > 0, "extent must be positive");
    }

    let dim = cell.dimension();
    let mut graph = Self::new(dim);
    let num_cells = extent.iter().fold(1usize, |acc, value| acc * (*value as usize));

    for cell_index in 0..num_cells {
      let cell_offset = index_to_offset(cell_index, extent);
      let cell_offset_f = CoordinateVector::from_iterator(
        dim,
        cell_offset.iter().map(|value| *value as f64),
      );
      for site in 0..cell.num_sites() {
        let coordinate = basis.basis_vectors() * (cell_offset_f.clone() + cell.site(site).coordinate.clone());
        graph.add_site(coordinate, cell.site(site).site_type);
      }
    }

    for cell_index in 0..num_cells {
      let cell_offset = index_to_offset(cell_index, extent);
      for bond_index in 0..cell.num_bonds() {
        let bond = cell.bond(bond_index);
        let mut target_offset = cell_offset.clone() + bond.target_offset.clone();
        if !wrap_offset(&mut target_offset, extent, boundary) {
          continue;
        }
        let target_cell = offset_to_index(&target_offset, extent);
        let source_site = cell_index * cell.num_sites() + bond.source;
        let target_site = target_cell * cell.num_sites() + bond.target;
        if source_site != target_site {
          graph.add_bond(source_site, target_site, bond.bond_type);
        }
      }
    }

    graph
  }

  pub fn from_basis_unitcell_length(
    basis: &Basis,
    cell: &Unitcell,
    length: usize,
    boundary: Boundary,
  ) -> Self {
    let extent = ExtentVector::from_element(cell.dimension(), length as i64);
    let boundary = vec![boundary; cell.dimension()];
    Self::from_basis_unitcell_extent(basis, cell, &extent, &boundary)
  }

  pub fn dimension(&self) -> usize {
    self.dim
  }

  pub fn num_sites(&self) -> usize {
    self.sites.len()
  }

  pub fn site_type(&self, site: usize) -> i32 {
    self.sites[site].site_type
  }

  pub fn coordinate(&self, site: usize) -> &CoordinateVector {
    &self.coordinates[site]
  }

  pub fn num_neighbors(&self, site: usize) -> usize {
    self.sites[site].neighbors.len()
  }

  pub fn neighbor(&self, site: usize, neighbor: usize) -> usize {
    self.sites[site].neighbors[neighbor]
  }

  pub fn neighbor_bond(&self, site: usize, neighbor: usize) -> usize {
    self.sites[site].neighbor_bonds[neighbor]
  }

  pub fn num_bonds(&self) -> usize {
    self.bonds.len()
  }

  pub fn bond_type(&self, bond: usize) -> i32 {
    self.bonds[bond].bond_type
  }

  pub fn source(&self, bond: usize) -> usize {
    self.bonds[bond].source
  }

  pub fn target(&self, bond: usize) -> usize {
    self.bonds[bond].target
  }

  pub fn edge_sites(&self, bond: usize) -> (usize, usize) {
    let bond = &self.bonds[bond];
    (bond.source, bond.target)
  }

  pub fn add_site(&mut self, coordinate: CoordinateVector, site_type: i32) -> usize {
    assert_eq!(coordinate.len(), self.dim, "dimension mismatch");
    let index = self.sites.len();
    self.sites.push(Site { site_type, neighbors: Vec::new(), neighbor_bonds: Vec::new() });
    self.coordinates.push(coordinate);
    index
  }

  pub fn add_bond(&mut self, source: usize, target: usize, bond_type: i32) -> usize {
    if source >= self.num_sites() || target >= self.num_sites() {
      panic!("site index out of range");
    }
    if source == target {
      panic!("self loop is not allowed");
    }
    let index = self.bonds.len();
    self.bonds.push(Bond { source, target, bond_type });
    self.sites[source].neighbors.push(target);
    self.sites[source].neighbor_bonds.push(index);
    self.sites[target].neighbors.push(source);
    self.sites[target].neighbor_bonds.push(index);
    index
  }
}

fn index_to_offset(index: usize, extent: &ExtentVector) -> OffsetVector {
  let dim = extent.len();
  let mut remainder = index;
  let mut offset = OffsetVector::zeros(dim);
  for axis in 0..dim {
    let size = extent[axis] as usize;
    offset[axis] = (remainder % size) as i64;
    remainder /= size;
  }
  offset
}

fn offset_to_index(offset: &OffsetVector, extent: &ExtentVector) -> usize {
  let mut index = 0usize;
  let mut stride = 1usize;
  for axis in 0..extent.len() {
    index += (offset[axis] as usize) * stride;
    stride *= extent[axis] as usize;
  }
  index
}

fn wrap_offset(offset: &mut OffsetVector, extent: &ExtentVector, boundary: &[Boundary]) -> bool {
  for axis in 0..offset.len() {
    let size = extent[axis];
    let value = offset[axis];
    match boundary[axis] {
      Boundary::Open => {
        if value < 0 || value >= size {
          return false;
        }
      }
      Boundary::Periodic => {
        let wrapped = ((value % size) + size) % size;
        offset[axis] = wrapped;
      }
    }
  }
  true
}
