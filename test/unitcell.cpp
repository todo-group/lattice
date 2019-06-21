#include "gtest/gtest.h"
#include "lattice/unitcell.hpp"

TEST(UnitcellTest, SimpleSquare) {
  lattice::basis basis("square lattice", 2);
  lattice::unitcell unitcell(basis);
  lattice::coordinate_t pos(2); pos << 0.0, 0.0;
  auto s0 = unitcell.add_vertex(pos, 0);
  lattice::offset_t offset_x(2); offset_x << 1, 0;
  auto b0 = unitcell.add_edge(s0, s0, offset_x, 0);
  lattice::offset_t offset_y(2); offset_y << 0, 1;
  auto b1 = unitcell.add_edge(s0, s0, offset_y, 0);

  EXPECT_EQ(2, unitcell.dimension());
  EXPECT_EQ(1, unitcell.num_vertices());
  EXPECT_EQ(2, unitcell.num_edges());
  EXPECT_EQ(0, unitcell.edge(0).source);
  EXPECT_EQ(0, unitcell.edge(0).target);
  EXPECT_EQ(0, unitcell.edge(1).source);
  EXPECT_EQ(0, unitcell.edge(1).target);
  // EXPECT_ANY_THROW(unitcell.add_edge(0, 1, offset_x, 0));
}
