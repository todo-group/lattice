import lattice


def describe(name, graph):
    print(f"{name}:")
    print(f"  dimension = {graph.dimension}")
    print(f"  sites     = {graph.num_sites}")
    print(f"  bonds     = {graph.num_bonds}")
    print(f"  first edge = {graph.edge_sites(0)}")


def main():
    chain = lattice.Graph.simple(1, 16)
    describe("periodic chain", chain)

    square = lattice.Graph.simple(2, 4)
    describe("periodic square lattice", square)
    print(f"  coordinates shape = {square.coordinates().shape}")
    print(f"  edges shape       = {square.edges().shape}")

    basis = lattice.Basis([[1.0, 0.0], [0.0, 1.0]])
    unitcell = lattice.Unitcell(2)
    unitcell.add_site([0.0, 0.0])
    unitcell.add_bond(0, 0, [1, 0])
    unitcell.add_bond(0, 0, [0, 1])

    generic_square = lattice.Graph.from_basis_unitcell_extent(
        basis,
        unitcell,
        [4, 4],
        [lattice.Boundary.Periodic, lattice.Boundary.Periodic],
    )
    describe("generic square lattice", generic_square)

    complete = lattice.Graph.fully_connected(10)
    describe("fully connected graph", complete)


if __name__ == "__main__":
    main()
