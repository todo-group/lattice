#ifndef LATTICE_RUST_XML_UNITCELL_HPP
#define LATTICE_RUST_XML_UNITCELL_HPP

#include "unitcell.hpp"
#include "rust_xml.hpp"

namespace lattice {

inline bool read_xml(ptree& pt, const std::string& name, unitcell& cell) {
  auto xml = detail::ptree_to_xml(pt);
  auto raw = lattice_unitcell_from_xml(xml.c_str(), name.c_str());
  if (!raw) {
    return false;
  }
  cell = unitcell(raw->dim);
  for (std::size_t s = 0; s < raw->num_sites; ++s) {
    coordinate_t pos = coordinate_t::Zero(raw->dim);
    for (std::size_t m = 0; m < raw->dim; ++m) pos(m) = raw->site_coordinates[s * raw->dim + m];
    cell.add_site(pos, raw->site_types[s]);
  }
  for (std::size_t b = 0; b < raw->num_bonds; ++b) {
    offset_t off = offset_t::Zero(raw->dim);
    for (std::size_t m = 0; m < raw->dim; ++m) off(m) = raw->bond_offsets[b * raw->dim + m];
    cell.add_bond(raw->bond_sources[b], raw->bond_targets[b], off, raw->bond_types[b]);
  }
  lattice_unitcell_raw_free(raw);
  return true;
}

inline ptree& write_xml(ptree& pt, const std::string& name, const unitcell& cell) {
  lattice_unitcell_raw raw{};
  raw.dim = cell.dimension();
  raw.num_sites = cell.num_sites();
  raw.num_bonds = cell.num_bonds();
  std::vector<int> site_types(raw.num_sites);
  std::vector<double> site_coordinates(raw.num_sites * raw.dim);
  for (std::size_t s = 0; s < cell.num_sites(); ++s) {
    site_types[s] = cell.site(s).type;
    for (std::size_t m = 0; m < raw.dim; ++m) site_coordinates[s * raw.dim + m] = cell.site(s).coordinate(m);
  }
  std::vector<std::size_t> bond_sources(raw.num_bonds);
  std::vector<std::size_t> bond_targets(raw.num_bonds);
  std::vector<int> bond_types(raw.num_bonds);
  std::vector<std::int64_t> bond_offsets(raw.num_bonds * raw.dim);
  for (std::size_t b = 0; b < cell.num_bonds(); ++b) {
    bond_sources[b] = cell.bond(b).source;
    bond_targets[b] = cell.bond(b).target;
    bond_types[b] = cell.bond(b).type;
    for (std::size_t m = 0; m < raw.dim; ++m) bond_offsets[b * raw.dim + m] = cell.bond(b).target_offset(m);
  }
  raw.site_types = site_types.data();
  raw.site_coordinates_len = site_coordinates.size();
  raw.site_coordinates = site_coordinates.data();
  raw.bond_sources = bond_sources.data();
  raw.bond_targets = bond_targets.data();
  raw.bond_types = bond_types.data();
  raw.bond_offsets_len = bond_offsets.size();
  raw.bond_offsets = bond_offsets.data();
  auto xml = detail::c_string(lattice_unitcell_to_xml(name.c_str(), &raw));
  ptree parsed;
  detail::xml_to_ptree(xml, parsed);
  pt.add_child("LATTICES.UNITCELL", parsed.get_child("LATTICES.UNITCELL"));
  return pt;
}

} // end namespace lattice

#endif
