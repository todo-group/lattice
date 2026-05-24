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

using namespace lattice;

class BasisIoTest : public testing::Test {
public:
  BasisIoTest() :
    basis1(basis::simple(1)),
    basis2(basis::simple(2)) {
  }
protected:
  basis basis1, basis2;
};

TEST_F(BasisIoTest, WriteXML) {
  const std::string xml = write_xml("simple1d", basis1);
  EXPECT_FALSE(xml.empty());
  EXPECT_NE(std::string::npos, xml.find("<LATTICES>"));
  EXPECT_NE(std::string::npos, xml.find("name=\"simple1d\""));
}

TEST_F(BasisIoTest, ReadXML) {
  const std::string xml = R"(
<LATTICES>
  <LATTICE name="simple1d" dimension="1">
    <BASIS><VECTOR>1</VECTOR></BASIS>
  </LATTICE>
  <LATTICE name="simple2d" dimension="2">
    <BASIS><VECTOR>1 0</VECTOR><VECTOR>0 1</VECTOR></BASIS>
  </LATTICE>
</LATTICES>
  )";

  basis bs;
  EXPECT_TRUE(read_xml(xml, "simple1d", bs));
  EXPECT_EQ(basis1.basis_vectors(), bs.basis_vectors());
  EXPECT_TRUE(read_xml(xml, "simple2d", bs));
  EXPECT_EQ(basis2.basis_vectors(), bs.basis_vectors());
}
