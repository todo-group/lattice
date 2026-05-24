use crate::types::{CoordinateVector, OffsetVector};

#[derive(Clone, Debug, PartialEq)]
pub struct Site {
  pub coordinate: CoordinateVector,
  pub site_type: i32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Bond {
  pub source: usize,
  pub target: usize,
  pub target_offset: OffsetVector,
  pub bond_type: i32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Unitcell {
  dim: usize,
  sites: Vec<Site>,
  bonds: Vec<Bond>,
}

impl Unitcell {
  pub fn new(dim: usize) -> Self {
    Self { dim, sites: Vec::new(), bonds: Vec::new() }
  }

  pub fn simple(dim: usize) -> Self {
    let mut cell = Self::new(dim);
    cell.add_site(CoordinateVector::zeros(dim), 0);
    for axis in 0..dim {
      let mut offset = OffsetVector::zeros(dim);
      offset[axis] = 1;
      cell.add_bond(0, 0, offset, 0);
    }
    cell
  }

  pub fn dimension(&self) -> usize {
    self.dim
  }

  pub fn num_sites(&self) -> usize {
    self.sites.len()
  }

  pub fn num_bonds(&self) -> usize {
    self.bonds.len()
  }

  pub fn site(&self, index: usize) -> &Site {
    &self.sites[index]
  }

  pub fn bond(&self, index: usize) -> &Bond {
    &self.bonds[index]
  }

  pub fn max_neighbors(&self) -> usize {
    let mut counts = vec![0usize; self.num_sites()];
    for bond in &self.bonds {
      counts[bond.source] += 1;
      counts[bond.target] += 1;
    }
    counts.into_iter().max().unwrap_or(0)
  }

  pub fn add_site(&mut self, coordinate: CoordinateVector, site_type: i32) -> usize {
    assert_eq!(coordinate.len(), self.dim, "site coordinate dimension mismatch");
    for value in coordinate.iter() {
      if *value < 0.0 || *value >= 1.0 {
        panic!("site coordinate out of range");
      }
    }
    let index = self.sites.len();
    self.sites.push(Site { coordinate, site_type });
    index
  }

  pub fn add_bond(&mut self, source: usize, target: usize, target_offset: OffsetVector, bond_type: i32) -> usize {
    if source >= self.num_sites() || target >= self.num_sites() {
      panic!("site index out of range");
    }
    if target_offset.len() != self.dim {
      panic!("unitcell offset dimension mismatch");
    }
    let index = self.bonds.len();
    self.bonds.push(Bond { source, target, target_offset, bond_type });
    index
  }
}
