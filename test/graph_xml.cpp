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

#include "gtest/gtest.h"
#include "lattice/graph_xml.hpp"

using namespace lattice;
using boost::property_tree::ptree;

TEST(GraphXMLTest, WriteXML) {
  ptree pt;
  ptree& root = pt.put("LATTICES", "");

  write_xml(root, "square lattice", graph::simple(2, 4));
  write_xml(std::cerr, pt,
    boost::property_tree::xml_writer_make_settings<std::string>(' ', 2));
}

TEST(GraphXMLTest, ReadXML1) {
  std::istringstream is(R"(
<GRAPH name="5-site dimerized" vertices="5">
  <EDGE type="0" source="1" target="2"/>
  <EDGE type="1" source="2" target="3"/>
  <EDGE type="0" source="3" target="4"/>
  <EDGE type="1" source="4" target="5"/>
</GRAPH>
    )");

  ptree pt;
  ptree& root = pt.put("LATTICES", "");
  read_xml(is, root);

  graph lat;
  EXPECT_TRUE(read_xml(pt, "5-site dimerized", lat));
  EXPECT_EQ(0, lat.dimension());
  EXPECT_EQ(5, lat.num_sites());
  EXPECT_EQ(4, lat.num_bonds());
  EXPECT_EQ(0, lat.bond_type(0));
  EXPECT_EQ(0, lat.source(0));
  EXPECT_EQ(1, lat.target(0));
  EXPECT_EQ(1, lat.bond_type(3));
  EXPECT_EQ(3, lat.source(3));
  EXPECT_EQ(4, lat.target(3));
}
