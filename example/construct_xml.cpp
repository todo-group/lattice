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

#include <iostream>
#include "lattice/graph_xml.hpp"

int main(int argc, char **argv) {
  std::string file = "lattices.xml";
  std::string basis_name = "square lattice";
  std::string cell_name = "simple2d";
  std::size_t length = 4;
  if (argc > 1) {
    if (argc == 4 || argc == 5) {
      file = argv[1];
      basis_name = argv[2];
      cell_name = argv[3];
      if (argc == 5) length = atoi(argv[4]);
    } else {
      std::cerr << "Error: " << argv[0] << " xmlfile basis cell [length]\n";
      exit(127);
    }
  }

  std::ifstream is(file);
  boost::property_tree::ptree pt;
  read_xml(is, pt);
  lattice::basis bs;
  read_xml(pt, basis_name, bs);
  lattice::unitcell cell;
  read_xml(pt, cell_name, cell);
  switch (cell.dimension()) {
  case 1:
    { lattice::graph lat(bs, cell, lattice::extent(length)); lat.print(std::cout); }
    break;
  case 2:
    { lattice::graph lat(bs, cell, lattice::extent(length, length)); lat.print(std::cout); }
    break;
  case 3:
    { lattice::graph lat(bs, cell, lattice::extent(length, length, length)); lat.print(std::cout); }
    break;
  default:
    std::cerr << "Unsupported lattice dimension\n";
    exit(127);
    break;
  }
}
