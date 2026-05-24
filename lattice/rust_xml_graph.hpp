#ifndef LATTICE_RUST_XML_GRAPH_HPP
#define LATTICE_RUST_XML_GRAPH_HPP

#include "graph.hpp"
#include "rust_xml.hpp"

namespace lattice {

inline bool read_xml(ptree& pt, const std::string& name, graph& lat) {
  auto xml = detail::ptree_to_xml(pt);
  auto raw = lattice_graph_from_xml(xml.c_str(), name.c_str());
  if (!raw) {
    return false;
  }
  lat = graph(raw->dim);
  for (std::size_t s = 0; s < raw->num_sites; ++s) {
    coordinate_t pos = coordinate_t::Zero(raw->dim);
    for (std::size_t m = 0; m < raw->dim; ++m) pos(m) = raw->site_coordinates[s * raw->dim + m];
    lat.add_site(pos, raw->site_types[s]);
  }
  for (std::size_t b = 0; b < raw->num_bonds; ++b) {
    lat.add_bond(raw->bond_sources[b], raw->bond_targets[b], raw->bond_types[b]);
  }
  lattice_graph_raw_free(raw);
  return true;
}

inline ptree& write_xml(ptree& pt, const std::string& name, const graph& lat) {
  lattice_graph_raw raw{};
  raw.dim = lat.dimension();
  raw.num_sites = lat.num_sites();
  raw.num_bonds = lat.num_bonds();
  std::vector<int> site_types(raw.num_sites);
  std::vector<double> site_coordinates(raw.num_sites * raw.dim);
  for (std::size_t s = 0; s < lat.num_sites(); ++s) {
    site_types[s] = lat.site_type(s);
    for (std::size_t m = 0; m < raw.dim; ++m) site_coordinates[s * raw.dim + m] = lat.coordinate(s)(m);
  }
  std::vector<std::size_t> bond_sources(raw.num_bonds);
  std::vector<std::size_t> bond_targets(raw.num_bonds);
  std::vector<int> bond_types(raw.num_bonds);
  for (std::size_t b = 0; b < lat.num_bonds(); ++b) {
    bond_sources[b] = lat.source(b);
    bond_targets[b] = lat.target(b);
    bond_types[b] = lat.bond_type(b);
  }
  raw.site_types = site_types.data();
  raw.site_coordinates_len = site_coordinates.size();
  raw.site_coordinates = site_coordinates.data();
  raw.bond_sources = bond_sources.data();
  raw.bond_targets = bond_targets.data();
  raw.bond_types = bond_types.data();
  auto xml = detail::c_string(lattice_graph_to_xml(name.c_str(), &raw));
  ptree parsed;
  detail::xml_to_ptree(xml, parsed);
  pt.add_child("LATTICES.GRAPH", parsed.get_child("LATTICES.GRAPH"));
  return pt;
}

} // end namespace lattice

#endif
