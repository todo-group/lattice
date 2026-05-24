use lattice_core::Graph;
use std::path::PathBuf;

pub fn print_graph(graph: &Graph) {
    println!("dimension = {}", graph.dimension());
    println!("sites = {}", graph.num_sites());
    println!("bonds = {}", graph.num_bonds());

    for site in 0..graph.num_sites() {
        let coordinate = graph.coordinate(site);
        let coordinate_text = if coordinate.is_empty() {
            "[]".to_string()
        } else {
            let joined = coordinate
                .iter()
                .map(|value| format!("{value:.6}"))
                .collect::<Vec<_>>()
                .join(", ");
            format!("[{joined}]")
        };

        let neighbors = (0..graph.num_neighbors(site))
            .map(|k| graph.neighbor(site, k).to_string())
            .collect::<Vec<_>>()
            .join(", ");

        println!(
            "site {site}: type={} coord={} neighbors=[{}]",
            graph.site_type(site),
            coordinate_text,
            neighbors
        );
    }

    for bond in 0..graph.num_bonds() {
        println!(
            "bond {bond}: {} <-> {} type={}",
            graph.source(bond),
            graph.target(bond),
            graph.bond_type(bond)
        );
    }
}

#[allow(dead_code)]
pub fn resolve_lattices_xml() -> Result<PathBuf, String> {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let repo_root = manifest_dir
        .parent()
        .and_then(|path| path.parent())
        .ok_or_else(|| "failed to locate repository root from CARGO_MANIFEST_DIR".to_string())?;

    let candidates = [
        PathBuf::from("cxx/example/lattices.xml"),
        PathBuf::from("example/lattices.xml"),
        PathBuf::from("lattices.xml"),
        repo_root.join("cxx/example/lattices.xml"),
    ];

    for candidate in candidates {
        if candidate.exists() {
            return Ok(candidate);
        }
    }

    Err("could not find lattices.xml (tried cxx/example/lattices.xml, example/lattices.xml, and lattices.xml)".to_string())
}
