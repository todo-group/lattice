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
#include "lattice/unitcell.hpp"

using namespace lattice;

TEST(UnitcellTest, SimpleSquare1) {
  unitcell unitcell(2);
  coordinate_t pos(2); pos << 0.0, 0.0;
  auto s0 = unitcell.add_site(pos, 0);
  offset_t offset_x(2); offset_x << 1, 0;
  unitcell.add_bond(s0, s0, offset_x, 0);
  offset_t offset_y(2); offset_y << 0, 1;
  unitcell.add_bond(s0, s0, offset_y, 0);

  EXPECT_EQ(2, unitcell.dimension());
  EXPECT_EQ(1, unitcell.num_sites());
  EXPECT_EQ(2, unitcell.num_bonds());
  EXPECT_EQ(0, unitcell.bond(0).source);
  EXPECT_EQ(0, unitcell.bond(0).target);
  EXPECT_EQ(0, unitcell.bond(1).source);
  EXPECT_EQ(0, unitcell.bond(1).target);
  EXPECT_EQ(4, unitcell.max_neighbors());
  EXPECT_ANY_THROW(unitcell.add_bond(0, 1, offset_x, 0));
}

TEST(UnitcellTest, SimpleSquare2) {
  unitcell unitcell = unitcell::simple(2);

  EXPECT_EQ(2, unitcell.dimension());
  EXPECT_EQ(1, unitcell.num_sites());
  EXPECT_EQ(2, unitcell.num_bonds());
  EXPECT_EQ(0, unitcell.bond(0).source);
  EXPECT_EQ(0, unitcell.bond(0).target);
  EXPECT_EQ(0, unitcell.bond(1).source);
  EXPECT_EQ(0, unitcell.bond(1).target);
  EXPECT_EQ(4, unitcell.max_neighbors());
}
