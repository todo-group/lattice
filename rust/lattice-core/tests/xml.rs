use lattice_core::{read_basis_from_str, read_graph_from_str, read_unitcell_from_str, write_basis_to_string, write_graph_to_string, write_unitcell_to_string, Basis, CoordinateVector, ExtentVector, Graph, Unitcell};
use nalgebra::DMatrix;

const LATICES_XML: &str = include_str!("../../../example/lattices.xml");

#[test]
fn read_basis_from_repository_sample() {
  let basis = read_basis_from_str(LATICES_XML, "square lattice").expect("basis should parse");
  assert_eq!(basis.dimension(), 2);
  assert_eq!(basis.basis_vectors()[(0, 0)], 1.0);
  assert_eq!(basis.basis_vectors()[(1, 1)], 1.0);
}

#[test]
fn read_unitcell_from_repository_sample() {
  let cell = read_unitcell_from_str(LATICES_XML, "simple2d").expect("unitcell should parse");
  assert_eq!(cell.dimension(), 2);
  assert_eq!(cell.num_sites(), 1);
  assert_eq!(cell.num_bonds(), 2);
}

#[test]
fn read_graph_from_repository_sample() {
  let graph = read_graph_from_str(LATICES_XML, "triangle").expect("graph should parse");
  assert_eq!(graph.num_sites(), 3);
  assert_eq!(graph.num_bonds(), 3);
}

#[test]
fn basis_unitcell_and_graph_round_trip_to_xml() {
  let basis = Basis::new(DMatrix::identity(2, 2));
  let basis_xml = write_basis_to_string("square", &basis);
  let basis_again = read_basis_from_str(&basis_xml, "square").expect("basis should round-trip");
  assert_eq!(basis_again.dimension(), 2);

  let mut cell = Unitcell::new(2);
  cell.add_site(CoordinateVector::from_vec(vec![0.0, 0.0]), 0);
  let mut offset = ExtentVector::zeros(2);
  offset[0] = 1;
  cell.add_bond(0, 0, offset, 0);
  let cell_xml = write_unitcell_to_string("simple", &cell);
  let cell_again = read_unitcell_from_str(&cell_xml, "simple").expect("unitcell should round-trip");
  assert_eq!(cell_again.num_sites(), 1);
  assert_eq!(cell_again.num_bonds(), 1);

  let graph = Graph::simple(1, 4);
  let graph_xml = write_graph_to_string("chain", &graph);
  let graph_again = read_graph_from_str(&graph_xml, "chain").expect("graph should round-trip");
  assert_eq!(graph_again.num_sites(), 4);
  assert_eq!(graph_again.num_bonds(), 4);
}
