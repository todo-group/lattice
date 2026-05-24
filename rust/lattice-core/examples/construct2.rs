use lattice_core::{Basis, BasisMatrix, Boundary, CoordinateVector, ExtentVector, Graph, OffsetVector, Unitcell};

#[path = "support/mod.rs"]
mod support;

fn main() {
    let basis = Basis::new(BasisMatrix::from_row_slice(1, 1, &[1.0]));

    let mut unitcell = Unitcell::new(1);
    unitcell.add_site(CoordinateVector::from_element(1, 0.0), 0);
    unitcell.add_bond(0, 0, OffsetVector::from_element(1, 1), 0);

    let extent = ExtentVector::from_element(1, 16);
    let boundary = vec![Boundary::Periodic; 1];
    let graph = Graph::from_basis_unitcell_extent(&basis, &unitcell, &extent, &boundary);
    support::print_graph(&graph);
}
