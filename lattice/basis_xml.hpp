/*
   Copyright (C) 2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/lexical_cast.hpp>
#include "basis.hpp"

namespace lattice {

using boost::property_tree::ptree;

ptree& operator>>(ptree& pt, basis& bs) {
  std::string name;
  if (auto str = pt.get_optional<std::string>("<xmlattr>.name")) {
    name = str.get();
  }
  std::vector<std::vector<double>> basis_v;
  for (auto& child : pt.get_child("BASIS")) {
    if (child.first == "VECTOR") {
      basis_v.push_back(std::vector<double>());
      std::istringstream is(child.second.data());
      std::string s;
      while (is >> s) basis_v.back().push_back(boost::lexical_cast<double>(s));
    }
  }
  std::size_t dim;
  if (auto str = pt.get_optional<std::string>("<xmlattr>.dimension")) {
    dim = boost::lexical_cast<std::size_t>(str.get());
  } else {
    dim = basis_v.size();
  }
  // check dimension
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
  bs = basis(basis_in, name);
  return pt;
}

ptree& operator<<(ptree& pt, const basis& bs) {
  ptree& child = pt.add("LATTICE", "");
  child.put("<xmlattr>.name", bs.name());
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
  
} // end namespace lattice

#endif
