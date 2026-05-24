use lattice_core::Graph;

#[path = "support/mod.rs"]
mod support;

fn main() {
    let graph = Graph::simple(2, 4);
    support::print_graph(&graph);
}
