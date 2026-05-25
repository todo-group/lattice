**Python Bindings Plan**

Use **PyO3 + maturin** and add a new workspace crate:

`rust/lattice-python`

This should wrap `lattice-core` directly, not the existing `lattice-ffi`. The C ABI is useful for C++ compatibility, but Python will be cleaner, safer, and more ergonomic if it binds the Rust structs through PyO3.

**Phase 1: Packaging Skeleton**

Add:

- `rust/lattice-python/Cargo.toml`
- `rust/lattice-python/src/lib.rs`
- root or Python-package-level `pyproject.toml`
- workspace entry in `Cargo.toml`

Package name suggestion:

```toml
name = "lattice"
```

Rust crate name:

```toml
name = "lattice-python"
```

Python module:

```python
import lattice
```

Use maturin:

```sh
maturin develop
pytest
```

**Phase 2: Minimal Public API**

Expose Python classes matching the Rust model:

```python
Basis
Unitcell
Graph
Boundary
```

Initial constructors:

```python
Basis.simple(dim)
Basis(matrix)

Unitcell(dim)
Unitcell.simple(dim)

Graph(dim)
Graph.simple(dim, length)
Graph.fully_connected(num_sites)
Graph.from_basis_unitcell_extent(basis, unitcell, extent, boundary)
Graph.from_basis_unitcell_length(basis, unitcell, length, boundary)
```

Initial methods/properties:

```python
graph.dimension
graph.num_sites
graph.num_bonds
graph.site_type(i)
graph.coordinate(i)
graph.num_neighbors(i)
graph.neighbor(i, k)
graph.neighbor_bond(i, k)
graph.bond_type(i)
graph.source(i)
graph.target(i)
graph.edge_sites(i)
graph.add_site(coordinate, site_type=0)
graph.add_bond(source, target, bond_type=0)
```

For Python ergonomics, also add bulk accessors early:

```python
graph.coordinates()
graph.site_types()
graph.edges()
graph.bond_types()
```

Those will matter more than one-at-a-time calls for Python performance.

**Phase 3: XML API**

Expose the existing XML helpers as module functions:

```python
read_basis_from_string(xml, name)
read_unitcell_from_string(xml, name)
read_graph_from_string(xml, name)

read_basis_from_file(path, name)
read_unitcell_from_file(path, name)
read_graph_from_file(path, name)

write_basis_to_string(name, basis)
write_unitcell_to_string(name, unitcell)
write_graph_to_string(name, graph)
```

Maybe use Python naming aliases too:

```python
Basis.from_xml(...)
Graph.to_xml(...)
```

but keep the module-level functions first so they map clearly to Rust.

**Phase 4: NumPy Support**

Make NumPy optional but strongly recommended.

Return:

- basis matrix as `numpy.ndarray`
- graph coordinates as shape `(num_sites, dim)`
- edges as shape `(num_bonds, 2)`
- site/bond types as integer arrays

Accept Python sequences first, then optimize for NumPy inputs once the API is stable.

**Phase 5: Errors**

Replace Rust panics crossing the Python boundary with Python exceptions.

Important because current Rust methods use `assert!` / `panic!` for things like dimension mismatch and invalid site indices. The Python wrapper should validate before calling, returning:

```python
ValueError
IndexError
RuntimeError
```

Longer-term, `lattice-core` itself could move from panics to `Result` for public constructors/mutators, but the Python binding can shield users initially.

**Phase 6: Tests And Docs**

Add Python tests under something like:

```text
python/tests/
```

or:

```text
rust/lattice-python/tests/
```

Test first:

- `Graph.simple(1, 16)`
- `Graph.simple(2, 4)`
- `Unitcell.simple`
- XML read/write round trips
- invalid dimension/index errors
- NumPy array shapes if NumPy is enabled

Add README section:

```sh
pip install maturin
maturin develop -m rust/lattice-python/Cargo.toml
python -c "import lattice; print(lattice.Graph.simple(2, 4).num_sites)"
```

**Phase 7: Wheels / CI**

Once local bindings work:

- add GitHub Actions job using `maturin-action`
- build wheels for macOS, Linux, Windows
- test import from built wheel
- optionally publish to PyPI later

**Recommended First Milestone**

Implement a thin, usable MVP:

```python
import lattice

g = lattice.Graph.simple(2, 4)

assert g.dimension == 2
assert g.num_sites == 16
assert g.num_bonds == 32
print(g.coordinate(0))
print(g.edges())
```

That gives you a real Python package quickly, exercises the Rust core, and sets up the shape for Julia bindings later without disturbing the C++ bridge.
