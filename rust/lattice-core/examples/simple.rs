use lattice_core::Graph;

fn main() {
    let graph = Graph::simple(2, 4);
    println!("dimension={}", graph.dimension());
    println!("sites={}", graph.num_sites());
    println!("bonds={}", graph.num_bonds());
}
