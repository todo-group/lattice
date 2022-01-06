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

#ifndef LATTICE_GRAPH_XML_HPP
#define LATTICE_GRAPH_XML_HPP

#include <sstream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "graph.hpp"
#include "unitcell_xml.hpp"
#include "basis_xml.hpp"

namespace lattice {

using boost::property_tree::ptree;

inline ptree& operator>>(ptree& pt, graph& lat) {
  std::size_t dim = 0;
  if (auto str = pt.get_optional<std::string>("<xmlattr>.dimension")) {
    dim = stoi(str.get());
  }
  lat = graph(dim);

  std::size_t ns = 0;
  if (auto str = pt.get_optional<std::string>("<xmlattr>.vertices")) ns = stoi(str.get());

  for (auto& v : pt) {
    if (v.first == "VERTEX") {
      bool found_pos = false;
      coordinate_t pos = coordinate_t::Zero(dim);
      int tp = 0;
      if (lat.dimension() > 0) {
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
      }
      if (auto str = v.second.get_optional<std::string>("<xmlattr>.type")) tp = stoi(str.get());
      lat.add_site(pos, tp);
    }
  }
     
  if (lat.num_sites() > 0) {
    if (ns > 0 && lat.num_sites() != ns)
      throw std::invalid_argument("inconsistent number of sites");
  } else {
    coordinate_t pos = coordinate_t::Zero(dim);
    for (std::size_t i = 0; i < ns; ++i) lat.add_site(pos, 0);
  }
  
  for (auto& e : pt) {
    if (e.first == "EDGE") {
      std::size_t source, target; // offset one
      int tp = 0;
      if (auto str = e.second.get_optional<std::string>("<xmlattr>.source")) {
        source = stoi(str.get());
      } else {
        throw std::invalid_argument("source attribute not found");
      }
      if (auto str = e.second.get_optional<std::string>("<xmlattr>.target")) {
        target = stoi(str.get());
      } else {
        throw std::invalid_argument("target attribute not found");
      }
      if (auto str = e.second.get_optional<std::string>("<xmlattr>.type")) tp = stoi(str.get());
      lat.add_bond(source - 1, target - 1, tp);
    }
  }
  return pt;
}

inline ptree& write_xml(ptree& pt, const std::string& name, const graph& lat) {
  ptree& root = pt.add("GRAPH", "");
  root.put("<xmlattr>.name", name);
  if (lat.dimension() > 0) root.put("<xmlattr>.dimension", lat.dimension());
  root.put("<xmlattr>.vertices", lat.num_sites());
  for (std::size_t s = 0; s < lat.num_sites(); ++s) {
    ptree& vertex = root.add("VERTEX", "");
    vertex.put("<xmlattr>.type", lat.site_type(s));
    if (lat.dimension() > 0) {
      std::ostringstream os;
      os << std::setprecision(std::numeric_limits<double>::max_digits10);
      auto pos = lat.coordinate(s);
      for (std::size_t m = 0; m < lat.dimension(); ++m) os << (m ? " " : "") << pos(m);
      vertex.add("COORDINATE", os.str());
    }
  }
  for (std::size_t b = 0; b < lat.num_bonds(); ++b) {
    ptree& edge = root.add("EDGE", "");
    edge.put("<xmlattr>.type", lat.bond_type(b));
    edge.put("<xmlattr>.source", lat.source(b) + 1);
    edge.put("<xmlattr>.target", lat.target(b) + 1);
  }
  return pt;
}

inline bool read_xml(ptree& pt, const std::string& name, graph& lat) {
  for (auto& child : pt.get_child("LATTICES")) {
    if (child.first == "GRAPH") {
      if (auto str = child.second.get_optional<std::string>("<xmlattr>.name")) {
        if (str.get() == name) {
          child.second >> lat;
          return true;
        }
      }
    }
  }
  return false;
}
  
} // end namespace lattice

#endif
