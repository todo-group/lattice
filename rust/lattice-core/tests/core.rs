use lattice_core::{Basis, Boundary, CoordinateVector, ExtentVector, Graph, OffsetVector, Unitcell};
use nalgebra::DMatrix;

#[test]
fn basis_and_unitcell_validate_dimensions() {
  let basis = Basis::new(DMatrix::identity(2, 2));
  assert_eq!(basis.dimension(), 2);
  assert_eq!(basis.volume(), 1.0);

  let mut cell = Unitcell::new(2);
  let site = CoordinateVector::from_vec(vec![0.0, 0.0]);
  assert_eq!(cell.add_site(site, 7), 0);
  let offset = OffsetVector::from_vec(vec![1, 0]);
  assert_eq!(cell.add_bond(0, 0, offset, 3), 0);
  assert_eq!(cell.num_sites(), 1);
  assert_eq!(cell.num_bonds(), 1);
}

#[test]
fn graph_fully_connected_has_expected_counts() {
  let graph = Graph::fully_connected(4);
  assert_eq!(graph.dimension(), 0);
  assert_eq!(graph.num_sites(), 4);
  assert_eq!(graph.num_bonds(), 6);
  assert_eq!(graph.num_neighbors(0), 3);
}

#[test]
fn graph_simple_builds_periodic_chain() {
  let graph = Graph::simple(1, 4);
  assert_eq!(graph.dimension(), 1);
  assert_eq!(graph.num_sites(), 4);
  assert_eq!(graph.num_bonds(), 4);
  assert_eq!(graph.site_type(0), 0);
  assert_eq!(graph.coordinate(0)[0], 0.0);
}

#[test]
fn graph_can_build_from_cell_and_extent() {
  let basis = Basis::new(DMatrix::identity(2, 2));
  let mut cell = Unitcell::new(2);
  cell.add_site(CoordinateVector::from_vec(vec![0.0, 0.0]), 1);
  let mut x = OffsetVector::zeros(2);
  x[0] = 1;
  cell.add_bond(0, 0, x, 2);
  let extent = ExtentVector::from_vec(vec![2, 3]);
  let graph = Graph::from_basis_unitcell_extent(&basis, &cell, &extent, &[Boundary::Periodic, Boundary::Periodic]);
  assert_eq!(graph.num_sites(), 6);
  assert_eq!(graph.num_bonds(), 6);
  assert_eq!(graph.site_type(0), 1);
}
