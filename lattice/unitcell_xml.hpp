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

#ifndef LATTICE_UNITCELL_XML_HPP
#define LATTICE_UNITCELL_XML_HPP

#include <sstream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "unitcell.hpp"

namespace lattice {

using boost::property_tree::ptree;

ptree& operator>>(ptree& pt, unitcell& cell) {
  std::size_t dim = 0;
  if (auto str = pt.get_optional<std::string>("<xmlattr>.dimension")) {
    dim = stoi(str.get());
  } else {
    throw std::invalid_argument("attribute dimension does not found");
  }
  cell = unitcell(dim);

  std::size_t ns = 0;
  if (auto str = pt.get_optional<std::string>("<xmlattr>.vertices")) ns = stoi(str.get());

  for (auto& v : pt) {
    if (v.first == "VERTEX") {
      bool found_pos = false;
      coordinate_t pos = coordinate_t::Zero(dim);
      int tp = 0;
      for (auto& c : v.second) {
        if (c.first == "COORDINATE") {
          if (found_pos) throw std::invalid_argument("duplicated <COORDINATE> tag");
          std::istringstream is(c.second.data());
          std::string sin;
          for (std::size_t m = 0; m < dim; ++m) {
            is >> sin;
            pos(m) = stod(sin);
          }
          found_pos = true;
        }
      }
      if (auto str = v.second.get_optional<std::string>("<xmlattr>.type")) {
        tp = stoi(str.get());
      }
      cell.add_site(pos, tp);
    }
  }
     
  if (cell.num_sites() > 0) {
    if (ns > 0 && cell.num_sites() != ns)
      throw std::invalid_argument("inconsistent number of sites");
  } else {
    coordinate_t pos = coordinate_t::Zero(dim);
    for (std::size_t i = 0; i < ns; ++i) cell.add_site(pos, 0);
  }
  
  for (auto& e : pt) {
    if (e.first == "EDGE") {
      bool found_source = false, found_target = false;
      std::size_t source = 1, target = 1; // offset one
      offset_t source_offset = offset_t::Zero(dim), target_offset = offset_t::Zero(dim);
      int tp = 0;
      for (auto& st : e.second) {
        if (st.first == "SOURCE") {
          if (found_source) throw std::invalid_argument("duplicated <SOURCE> tag");
          if (auto str = st.second.get_optional<std::string>("<xmlattr>.vertex"))
            source = stoi(str.get());
          if (auto str = st.second.get_optional<std::string>("<xmlattr>.offset")) {
            std::istringstream is(str.get());
            std::string sin;
            for (std::size_t m = 0; m < dim; ++m) {
              is >> sin;
              source_offset(m) = stod(sin);
            }
          }
          found_source = true;
        } else if (st.first == "TARGET") {
          if (found_target) throw std::invalid_argument("duplicated <TARGET> tag");
          if (auto str = st.second.get_optional<std::string>("<xmlattr>.vertex"))
            target = stoi(str.get());
          if (auto str = st.second.get_optional<std::string>("<xmlattr>.offset")) {
            std::istringstream is(str.get());
            std::string sin;
            for (std::size_t m = 0; m < dim; ++m) {
              is >> sin;
              target_offset(m) = stod(sin);
            }
          }
          found_target = true;
        }
      }
      target_offset = target_offset - source_offset;
      if (auto str = e.second.get_optional<std::string>("<xmlattr>.type")) tp = stoi(str.get());
      cell.add_bond(source - 1, target - 1, target_offset, tp);
    }
  }
  return pt;
}

ptree& write_xml(ptree& pt, const std::string& name, const unitcell& cell) {
  ptree& root = pt.add("UNITCELL", "");
  root.put("<xmlattr>.name", name);
  root.put("<xmlattr>.dimension", cell.dimension());
  root.put("<xmlattr>.vertices", cell.num_sites());
  for (std::size_t s = 0; s < cell.num_sites(); ++s) {
    ptree& vertex = root.add("VERTEX", "");
    vertex.put("<xmlattr>.type", cell.site(s).type);
    std::ostringstream os;
    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    auto pos = cell.site(s).coordinate;
    for (std::size_t m = 0; m < cell.dimension(); ++m) os << (m ? " " : "") << pos(m);
    vertex.add("COORDINATE", os.str());
  }
  for (std::size_t b = 0; b < cell.num_bonds(); ++b) {
    ptree& edge = root.add("EDGE", "");
    edge.put("<xmlattr>.type", cell.bond(b).type);
    ptree& source = edge.add("SOURCE", "");
    source.put("<xmlattr>.source", cell.bond(b).source + 1);
    ptree& target = edge.add("TARGET", "");
    target.put("<xmlattr>.target", cell.bond(b).target + 1);
    std::ostringstream os;
    auto offset = cell.bond(b).target_offset;
    for (std::size_t m = 0; m < cell.dimension(); ++m) os << (m ? " " : "") << offset(m);
    target.put("<xmlattr>.offset", os.str());
  }
  return pt;
}

bool read_xml(ptree& pt, const std::string& name, unitcell& cell) {
  for (auto& child : pt.get_child("LATTICES")) {
    if (child.first == "UNITCELL") {
      if (auto str = child.second.get_optional<std::string>("<xmlattr>.name")) {
        if (str.get() == name) {
          child.second >> cell;
          return true;
        }
      }
    }
  }
  return false;
}
  
} // end namespace lattice

#endif
