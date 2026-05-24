#ifndef LATTICE_RUST_XML_FFI_H
#define LATTICE_RUST_XML_FFI_H

#include <cstddef>
#include <cstdint>

extern "C" {

struct lattice_basis_raw {
  std::size_t dim;
  std::size_t values_len;
  double* values;
};

struct lattice_unitcell_raw {
  std::size_t dim;
  std::size_t num_sites;
  int* site_types;
  std::size_t site_coordinates_len;
  double* site_coordinates;
  std::size_t num_bonds;
  std::size_t* bond_sources;
  std::size_t* bond_targets;
  int* bond_types;
  std::size_t bond_offsets_len;
  std::int64_t* bond_offsets;
};

struct lattice_graph_raw {
  std::size_t dim;
  std::size_t num_sites;
  int* site_types;
  std::size_t site_coordinates_len;
  double* site_coordinates;
  std::size_t num_bonds;
  std::size_t* bond_sources;
  std::size_t* bond_targets;
  int* bond_types;
};

char* lattice_basis_to_xml(const char* name, const lattice_basis_raw* raw);
lattice_basis_raw* lattice_basis_from_xml(const char* xml, const char* name);
void lattice_basis_raw_free(lattice_basis_raw* raw);

char* lattice_unitcell_to_xml(const char* name, const lattice_unitcell_raw* raw);
lattice_unitcell_raw* lattice_unitcell_from_xml(const char* xml, const char* name);
void lattice_unitcell_raw_free(lattice_unitcell_raw* raw);

char* lattice_graph_to_xml(const char* name, const lattice_graph_raw* raw);
lattice_graph_raw* lattice_graph_from_xml(const char* xml, const char* name);
void lattice_graph_raw_free(lattice_graph_raw* raw);

void lattice_string_free(char* ptr);

}

#endif
