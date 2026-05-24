#ifndef LATTICE_RUST_XML_HPP
#define LATTICE_RUST_XML_HPP

#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "basis.hpp"
#include "rust_xml_ffi.h"

namespace lattice {

using boost::property_tree::ptree;

namespace detail {

inline std::string ptree_to_xml(const ptree& pt) {
  std::ostringstream os;
  boost::property_tree::write_xml(os, pt);
  return os.str();
}

inline void xml_to_ptree(const std::string& xml, ptree& pt) {
  std::istringstream is(xml);
  boost::property_tree::read_xml(is, pt);
}

inline std::string c_string(char* ptr) {
  std::string value = ptr ? ptr : "";
  lattice_string_free(ptr);
  return value;
}

} // namespace detail

inline bool read_xml(ptree& pt, const std::string& name, basis& bs) {
  auto xml = detail::ptree_to_xml(pt);
  auto raw = lattice_basis_from_xml(xml.c_str(), name.c_str());
  if (!raw) {
    return false;
  }
  basis_t basis_in(raw->dim, raw->dim);
  for (std::size_t row = 0; row < raw->dim; ++row) {
    for (std::size_t col = 0; col < raw->dim; ++col) {
      basis_in(row, col) = raw->values[row * raw->dim + col];
    }
  }
  bs = basis(basis_in);
  lattice_basis_raw_free(raw);
  return true;
}

inline ptree& write_xml(ptree& pt, const std::string& name, const basis& bs) {
  lattice_basis_raw raw{};
  raw.dim = bs.dimension();
  raw.values_len = raw.dim * raw.dim;
  std::vector<double> values(raw.values_len);
  for (std::size_t row = 0; row < raw.dim; ++row) {
    for (std::size_t col = 0; col < raw.dim; ++col) {
      values[row * raw.dim + col] = bs.basis_vectors()(row, col);
    }
  }
  raw.values = values.data();
  auto xml = detail::c_string(lattice_basis_to_xml(name.c_str(), &raw));
  ptree parsed;
  detail::xml_to_ptree(xml, parsed);
  pt.add_child("LATTICES.LATTICE", parsed.get_child("LATTICES.LATTICE"));
  return pt;
}

} // end namespace lattice

#endif
