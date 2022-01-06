/*
   Copyright (C) 2019-2022 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#ifndef LATTICE_BASIS_XML_HPP
#define LATTICE_BASIS_XML_HPP

#include <sstream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "basis.hpp"

namespace lattice {

using boost::property_tree::ptree;

inline ptree& operator>>(ptree& pt, basis& bs) {
  std::vector<std::vector<double>> basis_v;
  for (auto& child : pt.get_child("BASIS")) {
    if (child.first == "VECTOR") {
      basis_v.push_back(std::vector<double>());
      std::istringstream is(child.second.data());
      std::string s;
      while (is >> s) basis_v.back().push_back(stod(s));
    }
  }
  std::size_t dim;
  if (auto str = pt.get_optional<std::string>("<xmlattr>.dimension")) {
    dim = stoi(str.get());
  } else {
    dim = basis_v.size();
  }
  if (dim == 0 || basis_v.size() != dim) {
    throw std::runtime_error("basis dimension mismatch");
  }
  for (std::size_t i = 0; i < dim; ++i) {
    if (basis_v[i].size() != dim) {
      throw std::runtime_error("basis dimension mismatch");
    }
  }
  basis_t basis_in(dim, dim);
  for (std::size_t i = 0; i < dim; ++i) {
    for (std::size_t j = 0; j < dim; ++j) {
      basis_in(j, i) = basis_v[i][j];
    }
  }
  bs = basis(basis_in);
  return pt;
}

inline ptree& write_xml(ptree& pt, const std::string& name, const basis& bs) {
  ptree& child = pt.add("LATTICE", "");
  child.put("<xmlattr>.name", name);
  child.put("<xmlattr>.dimension", bs.dimension());
  ptree& basis = child.add("BASIS", "");
  for (std::size_t j = 0; j < bs.dimension(); ++j) {
    std::ostringstream os;
    for (std::size_t i = 0; i < bs.dimension(); ++i) {
      os << (i ? " " : "") << bs.basis_vectors()(i, j);
    }
    basis.add("VECTOR", os.str());
  }
  return pt;
}

inline bool read_xml(ptree& pt, const std::string& name, basis& bs) {
  for (auto& child : pt.get_child("LATTICES")) {
    if (child.first == "LATTICE") {
      if (auto str = child.second.get_optional<std::string>("<xmlattr>.name")) {
        if (str.get() == name) {
          child.second >> bs;
          return true;
        }
      }
    }
  }
  return false;
}
  
} // end namespace lattice

#endif
