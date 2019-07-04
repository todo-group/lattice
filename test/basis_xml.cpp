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
#include "lattice/basis_xml.hpp"

using boost::property_tree::ptree;

class BasisIoTest : public testing::Test {
public:
  BasisIoTest() : basis1(1, "chain lattice"), basis2(2, "square lattice") {}
protected:
  lattice::basis basis1, basis2;
};

TEST_F(BasisIoTest, WriteXML) {
  ptree pt;
  ptree& root = pt.put("LATTICES", "");

  root << basis1 << basis2;
  write_xml(std::cerr, pt,
    boost::property_tree::xml_writer_make_settings<std::string>(' ', 2));
}

TEST_F(BasisIoTest, ReadXML) {
  ptree pt;
  ptree& root = pt.put("LATTICES", "");

  root << basis1 << basis2;

  std::map<std::string, lattice::basis> lattices;
  for (auto& child : pt.get_child("LATTICES")) {
    if (child.first == "LATTICE") {
      lattice::basis bs;
      child.second >> bs;
      lattices[bs.name()] = bs;
    }
  }
  EXPECT_EQ(basis1.basis_vectors(), lattices[basis1.name()].basis_vectors());
  EXPECT_EQ(basis2.basis_vectors(), lattices[basis2.name()].basis_vectors());
}
