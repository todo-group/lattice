use lattice_core::Graph;

#[path = "support/mod.rs"]
mod support;

fn main() {
    let graph = Graph::simple(1, 16);
    support::print_graph(&graph);
}
