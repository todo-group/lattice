![](https://github.com/todo-group/lattice/workflows/build/badge.svg)

# lattice

Simple Lattice/Graph Library

## Main functions

* Access to information of lattice structure (sites, bonds, etc)
* Construct lattice from unitcell, span vector of supercell, and boundary conditions
* Construct lattice by adding sites and bonds one by one
* Reading and writing ALPS Lattice XML file

## Prerequisites

* C++-14 compiler
* CMake to build tests and examples
* Eigen3
* [optional] Boost Library (Property Tree) for reading and writing XML files

## Rust workspace (work in progress)

The repository is being extended with a Rust core under `rust/` as the shared implementation base for future Python and Julia bindings.

### C++ XML compatibility bridge (opt-in)

You can opt-in to the Rust-backed XML implementation from C++ headers by enabling:

* `-DLATTICE_USE_RUST_XML=ON`

To also build the Rust C ABI library from CMake, enable:

* `-DLATTICE_BUILD_RUST_FFI=ON`

Rust targets are managed in the workspace under `rust/`:

* `lattice-core`: core model + XML parser/writer
* `lattice-ffi`: C ABI layer for C++ compatibility

## Build, test, and sample run

### C++ (default)

Configure and build:

```
mkdir build
cd build
cmake ..
cmake --build .
```

Run C++ tests:

```
ctest --output-on-failure
```

Run C++ samples:

```
./example/construct1
./example/construct2
./example/construct3
./example/construct4
./example/construct_xml
./example/ising
```

### Recommended on macOS (fixed SDK)

Use CMake presets to pin the SDK to `MacOSX15.4.sdk`:

```
cmake --preset macos-sdk154
cmake --build --preset macos-sdk154
ctest --preset macos-sdk154
```

### Rust XML smoke test (C++ compatibility path)

```
cmake --preset macos-sdk154-rustxml
cmake --build --preset macos-sdk154-rustxml-smoke
./build-sdk154-rustxml/test/rust_xml_bridge
```

### Rust

Build and run Rust tests:

```
cd rust
cargo build
cargo test
```

Run Rust sample:

```
cd rust
cargo run -p lattice-core --example simple
```

## Classes/types

* lattice::basis

  Helper class that contains the shape, i.e. the set of basis vectors, of the unit cell.

* lattice::unitcell

  Helper class that contains the structure, i.e. sites and bonds, of the unit cell.
  
* lattice::graph

  This class contains the structure of the whole lattice structure. It provides various information of sites (vertices) and bonds (edges) via the following member functions:

  |  member functions  |  description  |
  | ---- | ---- |
  | std::size\_t dimension() const; | dimension of the lattice |
  |||
  | std::size\_t num\_sites() const; | total number of sites |
  | std::size\_t site\_type(std::size\_t s) const; | type of site s |
  | const coordinate\_t& coordinate(std::size\_t s) const; | coordinate of site s |
  | std::size\_t num\_neighbors(std::size\_t s) const; | number of neighboring sites of site s  |
  | std::size\_t neighbor(std::size\_t s, std::size\_t k) const; | k-th neighbor site of site s |
  | std::size\_t neighbor\_bond(std::size\_t s, std::size\_t k) const; | bond connecting site s and its k-th neighbor site |
  |||
  | std::size\_t num\_bonds() const; | total number of bonds |
  | int bond\_type(std::size\_t b) const; | type of bond s |
  | std::size\_t source(std::size\_t b) const; | start point (site) of bond b |
  | std::size\_t target(std::size\_t b) const; | end point (site) of bond b |
  
## How to construct lattices

* periodic chain lattice of 16 sites

  * simplest interface

     ```
     lattice::graph lat = lattice::graph::simple(1, 16);
     ```
    
  * most generic interface

     ```
     lattice::basis_t bs(1, 1); bs << 1; // 1x1 matrix
     lattice::basis basis(bs);
     lattice::unitcell unitcell(1);
     unitcell.add_site(lattice::coordinate(0), 0);
     unitcell.add_bond(0, 0, lattice::offset(1), 0);
     lattice::span_t span(1, 1); span << 16; // 1x1 matrix
     std::vector<lattice::boundary_t> boundary(1, lattice::boundary_t::periodic);
     lattice::graph lat(basis, unitcell, span, boundary);
     ```
     
* periodic square lattice of 4 x 4 sites

  * simplest interface

     ```
     lattice::graph lat = lattice::graph::simple(2, 4);
     ```
    
  * most generic interface

     ```
     lattice::basis_t bs(2, 2); bs << 1, 0, 0, 1; // 2x2 matrix
     lattice::basis basis(bs);
     lattice::unitcell unitcell(2);
     unitcell.add_site(lattice::coordinate(0, 0), 0);
     unitcell.add_bond(0, 0, lattice::offset(1, 0), 0);
     unitcell.add_bond(0, 0, lattice::offset(0, 1), 0);
     lattice::span_t span(2, 2); span << 4, 0, 0, 4; // 2x2 matrix
     std::vector<lattice::boundary_t> boundary(2, lattice::boundary_t::periodic);
     lattice::graph lat(basis, unitcell, span, boundary);
     ```

   * reading basis and unitcell from XML file
   
      ```
      std::string file = "lattices.xml";
      std::ifstream is(file);
      boost::property_tree::ptree pt;
      read_xml(is, pt);
      lattice::basis bs;
      read_xml(pt, "square lattice", bs);
      lattice::unitcell cell;
      read_xml(pt, "simple2d", cell);
      lattice::graph lat(bs, cell, lattice::extent(4, 4));
      ```

* fully connected lattice of 10 sites

  ```
  lattice::graph lat = lattice::graph::fully_connected(10);
  ```
