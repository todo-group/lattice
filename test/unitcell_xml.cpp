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

#include <sstream>
#include "gtest/gtest.h"
#include "lattice/unitcell_xml.hpp"

using namespace lattice;
using boost::property_tree::ptree;

TEST(UnitcellXMLTest, WriteXML) {
  ptree pt;
  ptree& root = pt.put("LATTICES", "");

  write_xml(root, "simple1d", unitcell::simple(1));
  write_xml(root, "simple2d", unitcell::simple(2));
  write_xml(std::cerr, pt,
    boost::property_tree::xml_writer_make_settings<std::string>(' ', 2));
  
  unitcell cell1, cell2;
  EXPECT_TRUE(read_xml(pt, "simple1d", cell1));
  EXPECT_TRUE(read_xml(pt, "simple1d", cell2));
}

TEST(UnitcellXMLTest, ReadXML1) {
  std::istringstream is(R"(
<UNITCELL name="simple1d" dimension="1">
  <VERTEX/>
  <EDGE><SOURCE vertex="1"/><TARGET vertex="1" offset="1"/></EDGE>
</UNITCELL>
    )");

  ptree pt;
  ptree& root = pt.put("LATTICES", "");
  read_xml(is, root);

  unitcell cell;
  EXPECT_TRUE(read_xml(pt, "simple1d", cell));
  EXPECT_EQ(1, cell.dimension());
  EXPECT_EQ(1, cell.num_sites());
  EXPECT_EQ(1, cell.num_bonds());
  EXPECT_EQ(0, cell.bond(0).source);
  EXPECT_EQ(0, cell.bond(0).target);
  offset_t offset(1);
  offset << 1.0;
  EXPECT_EQ(offset, cell.bond(0).target_offset);
}

TEST(UnitcellXMLTest, ReadXML2) {
  std::istringstream is(R"(
<UNITCELL name="kagome" dimension="2">
  <VERTEX><COORDINATE>0   0</COORDINATE></VERTEX>
  <VERTEX><COORDINATE>0.5 0</COORDINATE></VERTEX>
  <VERTEX><COORDINATE>0 0.5</COORDINATE></VERTEX>
  <EDGE><SOURCE vertex="1"/><TARGET vertex="2"/></EDGE>
  <EDGE><SOURCE vertex="1"/><TARGET vertex="3"/></EDGE>
  <EDGE><SOURCE vertex="2"/><TARGET vertex="3"/></EDGE>
  <EDGE><SOURCE vertex="1"/><TARGET vertex="2" offset="-1 0"/></EDGE>
  <EDGE><SOURCE vertex="2"/><TARGET vertex="3" offset="1 -1"/></EDGE>
  <EDGE><SOURCE vertex="1"/><TARGET vertex="3" offset="0 -1"/></EDGE>
</UNITCELL>
    )");

  ptree pt;
  ptree& root = pt.put("LATTICES", "");
  read_xml(is, root);

  unitcell cell;
  EXPECT_TRUE(read_xml(pt, "kagome", cell));
  EXPECT_EQ(2, cell.dimension());
  EXPECT_EQ(3, cell.num_sites());
  EXPECT_EQ(6, cell.num_bonds());
  EXPECT_EQ(0, cell.bond(0).source);
  EXPECT_EQ(1, cell.bond(0).target);
  EXPECT_EQ(0, cell.bond(5).source);
  EXPECT_EQ(2, cell.bond(5).target);
  offset_t offset(2);
  offset << 0, -1;
  EXPECT_EQ(offset, cell.bond(5).target_offset);
}

TEST(UnitcellXMLTest, ReadXML3) {
  std::istringstream is(R"(
<UNITCELL name="anisotropic3d" dimension="3">
  <VERTEX/>
  <EDGE type="0"><SOURCE vertex="1"/><TARGET vertex="1" offset="1 0 0"/></EDGE>
  <EDGE type="1"><SOURCE vertex="1"/><TARGET vertex="1" offset="0 1 0"/></EDGE>
  <EDGE type="2"><SOURCE vertex="1"/><TARGET vertex="1" offset="0 0 1"/></EDGE>
</UNITCELL>
    )");

  ptree pt;
  ptree& root = pt.put("LATTICES", "");
  read_xml(is, root);

  unitcell cell;
  EXPECT_TRUE(read_xml(pt, "anisotropic3d", cell));
  EXPECT_EQ(3, cell.dimension());
  EXPECT_EQ(1, cell.num_sites());
  EXPECT_EQ(3, cell.num_bonds());
  EXPECT_EQ(0, cell.bond(0).source);
  EXPECT_EQ(0, cell.bond(0).target);
  EXPECT_EQ(0, cell.bond(2).source);
  EXPECT_EQ(0, cell.bond(2).target);
  EXPECT_EQ(2, cell.bond(2).type);
  offset_t offset(3);
  offset << 0, 1, 0;
  EXPECT_EQ(offset, cell.bond(1).target_offset);
}
