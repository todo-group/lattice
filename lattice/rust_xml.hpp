#ifndef LATTICE_RUST_XML_HPP
#define LATTICE_RUST_XML_HPP

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "basis.hpp"
#include "rust_xml_ffi.h"

namespace lattice {

namespace detail {

inline std::string c_string(char* ptr) {
  std::string value = ptr ? ptr : "";
  lattice_string_free(ptr);
  return value;
}

inline std::string read_stream(std::istream& is) {
  std::ostringstream os;
  os << is.rdbuf();
  return os.str();
}

inline std::string last_error_message() {
  return c_string(lattice_last_error_message());
}

} // namespace detail

inline bool read_xml(const std::string& xml, const std::string& name, basis& bs) {
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

inline bool read_xml(std::istream& is, const std::string& name, basis& bs) {
  return read_xml(detail::read_stream(is), name, bs);
}

inline bool read_xml_file(const std::string& path, const std::string& name, basis& bs) {
  std::ifstream is(path);
  if (!is) {
    return false;
  }
  return read_xml(is, name, bs);
}

inline std::string write_xml(const std::string& name, const basis& bs) {
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
  return detail::c_string(lattice_basis_to_xml(name.c_str(), &raw));
}

} // end namespace lattice

#endif
