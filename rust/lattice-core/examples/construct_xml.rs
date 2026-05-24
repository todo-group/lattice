use lattice_core::{read_basis_from_file, read_unitcell_from_file, Boundary, ExtentVector, Graph};

#[path = "support/mod.rs"]
mod support;

fn main() {
    let args = std::env::args().collect::<Vec<_>>();

    let mut file = support::resolve_lattices_xml().unwrap_or_else(|message| {
        eprintln!("Error: {message}");
        std::process::exit(127);
    });
    let mut basis_name = "square lattice".to_string();
    let mut cell_name = "simple2d".to_string();
    let mut length: usize = 4;

    if args.len() > 1 {
        if args.len() == 4 || args.len() == 5 {
            file = std::path::PathBuf::from(&args[1]);
            basis_name = args[2].clone();
            cell_name = args[3].clone();
            if args.len() == 5 {
                length = args[4].parse::<usize>().unwrap_or_else(|_| {
                    eprintln!("Error: invalid length: {}", args[4]);
                    std::process::exit(127);
                });
            }
        } else {
            eprintln!("Error: {} xmlfile basis cell [length]", args[0]);
            std::process::exit(127);
        }
    }

    let basis = read_basis_from_file(&file, &basis_name).unwrap_or_else(|error| {
        eprintln!("Failed to read basis XML entry '{basis_name}': {error}");
        std::process::exit(127);
    });

    let cell = read_unitcell_from_file(&file, &cell_name).unwrap_or_else(|error| {
        eprintln!("Failed to read unitcell XML entry '{cell_name}': {error}");
        std::process::exit(127);
    });

    let extent = ExtentVector::from_element(cell.dimension(), length as i64);
    let boundary = vec![Boundary::Periodic; cell.dimension()];
    let graph = Graph::from_basis_unitcell_extent(&basis, &cell, &extent, &boundary);

    support::print_graph(&graph);
}
