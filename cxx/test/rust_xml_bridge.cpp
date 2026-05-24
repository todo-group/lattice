/*
   Copyright (C) 2026 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

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
#include <sstream>
#include <stdexcept>
#include <string>
#include "lattice/rust_xml_ffi.h"

namespace {

void Check(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

void TestReadFromLatticesDocument() {
  const char* xml = R"(
<LATTICES>
  <LATTICE name="square lattice" dimension="2">
    <BASIS>
      <VECTOR>1 0</VECTOR>
      <VECTOR>0 1</VECTOR>
    </BASIS>
  </LATTICE>
  <UNITCELL name="simple2d" dimension="2">
    <VERTEX/>
    <EDGE><SOURCE vertex="1"/><TARGET vertex="1" offset="1 0"/></EDGE>
    <EDGE><SOURCE vertex="1"/><TARGET vertex="1" offset="0 1"/></EDGE>
  </UNITCELL>
  <GRAPH name="triangle" vertices="3">
    <EDGE type="0" source="1" target="2"/>
    <EDGE type="0" source="2" target="3"/>
    <EDGE type="0" source="3" target="1"/>
  </GRAPH>
</LATTICES>
  )";

  lattice_basis_raw* bs = lattice_basis_from_xml(xml, "square lattice");
  lattice_unitcell_raw* cell = lattice_unitcell_from_xml(xml, "simple2d");
  lattice_graph_raw* lat = lattice_graph_from_xml(xml, "triangle");

  Check(bs != nullptr, "failed to read basis");
  Check(cell != nullptr, "failed to read unitcell");
  Check(lat != nullptr, "failed to read graph");

  Check(bs->dim == 2, "basis dimension mismatch");
  Check(cell->dim == 2, "unitcell dimension mismatch");
  Check(cell->num_sites == 1, "unitcell site count mismatch");
  Check(cell->num_bonds == 2, "unitcell bond count mismatch");
  Check(lat->dim == 0, "graph dimension mismatch");
  Check(lat->num_sites == 3, "graph site count mismatch");
  Check(lat->num_bonds == 3, "graph bond count mismatch");

  lattice_basis_raw_free(bs);
  lattice_unitcell_raw_free(cell);
  lattice_graph_raw_free(lat);
}

void TestWriteThenReadGraph() {
  lattice_graph_raw raw{};
  raw.dim = 1;
  raw.num_sites = 4;
  raw.num_bonds = 4;

  int site_types[] = {0, 0, 0, 0};
  double site_coordinates[] = {0.0, 1.0, 2.0, 3.0};
  std::size_t bond_sources[] = {0, 1, 2, 3};
  std::size_t bond_targets[] = {1, 2, 3, 0};
  int bond_types[] = {0, 0, 0, 0};

  raw.site_types = site_types;
  raw.site_coordinates_len = 4;
  raw.site_coordinates = site_coordinates;
  raw.bond_sources = bond_sources;
  raw.bond_targets = bond_targets;
  raw.bond_types = bond_types;

  char* xml = lattice_graph_to_xml("chain4", &raw);
  Check(xml != nullptr, "failed to write graph XML");

  lattice_graph_raw* loaded = lattice_graph_from_xml(xml, "chain4");
  lattice_string_free(xml);

  Check(loaded != nullptr, "failed to read back graph");
  Check(loaded->dim == 1, "graph dimension after round-trip mismatch");
  Check(loaded->num_sites == 4, "graph site count after round-trip mismatch");
  Check(loaded->num_bonds == 4, "graph bond count after round-trip mismatch");
  lattice_graph_raw_free(loaded);
}

} // namespace

int main() {
  try {
    TestReadFromLatticesDocument();
    TestWriteThenReadGraph();
    std::cout << "rust_xml_bridge: OK" << std::endl;
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "rust_xml_bridge: FAILED: " << e.what() << std::endl;
    return 1;
  }
}
