/*****************************************************************************
*
* Copyright (C) 2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include "gtest/gtest.h"
#include "lattice/unitcell.hpp"

TEST(UnitcellTest, SimpleSquare) {
  lattice::basis basis(2);
  lattice::unitcell unitcell(basis);
  lattice::coordinate_t pos(2); pos << 0.0, 0.0;
  auto s0 = unitcell.add_site(pos, 0);
  lattice::offset_t offset_x(2); offset_x << 1, 0;
  auto b0 = unitcell.add_bond(s0, s0, offset_x, 0);
  lattice::offset_t offset_y(2); offset_y << 0, 1;
  auto b1 = unitcell.add_bond(s0, s0, offset_y, 0);

  EXPECT_EQ(2, unitcell.dimension());
  EXPECT_EQ(1, unitcell.num_sites());
  EXPECT_EQ(2, unitcell.num_bonds());
  EXPECT_EQ(0, unitcell.bond(0).source);
  EXPECT_EQ(0, unitcell.bond(0).target);
  EXPECT_EQ(0, unitcell.bond(1).source);
  EXPECT_EQ(0, unitcell.bond(1).target);
  EXPECT_DOUBLE_EQ(1.0, unitcell.volume());
  EXPECT_EQ(4, unitcell.max_neighbors());
  // EXPECT_ANY_THROW(unitcell.add_bond(0, 1, offset_x, 0));
}
