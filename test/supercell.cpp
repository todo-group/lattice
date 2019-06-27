#include "gtest/gtest.h"
#include "lattice/supercell.hpp"

TEST(SupercellTest, TwoD) {
  lattice::span_t span(2, 2);
  span << 4, 2,
          1, 3;
  lattice::supercell supercell(span);
  lattice::offset_t offset(2);
  
  // within_supercell
  offset << 0, 0; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 1, 1; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 2, 1; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 1; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 2, 2; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 2; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 4, 2; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 3; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 5, 3; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 1, 0; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 4, 1; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 5, 2; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 6, 3; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 0, -1; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 0; EXPECT_FALSE(supercell.within_supercell(offset));

  for (std::size_t i = 0; i < supercell.num_cells(); ++i) {
    EXPECT_EQ(i, supercell.lcord2index(supercell.index2lcord(i)));
  }
  for (std::size_t i = 0; i < supercell.num_cells(); ++i) {
    auto offset = supercell.lcord2offset(supercell.index2lcord(i));
    EXPECT_EQ(i, supercell.lcord2index(supercell.offset2lcord(offset)));
  }
  
  std::size_t source, target, target_r;
  lattice::offset_t crossing(2), crossing_r(2);
  source = 0; offset << 1, 0; target = 7; crossing << 0, -1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 3; offset << -1, -1; target = 8; crossing << 0, -1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 1; offset << 4, 0; target = 7; crossing << 1, -1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 7; offset << -2, 0; target = 3; crossing << -1, 1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
}

TEST(SupercellTest, ThreeD) {
  lattice::span_t span(3, 3);
  span << 4, 0, 0,
          2, 4, 0,
          0, 0, 4;
  lattice::supercell supercell(span);
  lattice::offset_t offset(3);

  for (std::size_t i = 0; i < supercell.num_cells(); ++i) {
    EXPECT_EQ(i, supercell.lcord2index(supercell.index2lcord(i)));
  }
  for (std::size_t i = 0; i < supercell.num_cells(); ++i) {
    auto offset = supercell.lcord2offset(supercell.index2lcord(i));
    EXPECT_EQ(i, supercell.lcord2index(supercell.offset2lcord(offset)));
  }
  
  std::size_t source, target, target_r;
  source = 0; offset << 1, 0, 0; target = 12;
  EXPECT_EQ(target, supercell.add_offset(source, offset).first);
  source = 0; offset << 0, 1, 0; target = 1;
  EXPECT_EQ(target, supercell.add_offset(source, offset).first);
  source = 0; offset << 0, 0, 1; target = 16;
  EXPECT_EQ(target, supercell.add_offset(source, offset).first);
  source = 0; offset << -1, 0, 0; target = 7;
  EXPECT_EQ(target, supercell.add_offset(source, offset).first);
  source = 0; offset << 0, -1, 0; target = 8;
  EXPECT_EQ(target, supercell.add_offset(source, offset).first);
  source = 0; offset << 0, 0, -1; target = 48;
  EXPECT_EQ(target, supercell.add_offset(source, offset).first);
  source = 14; offset << 1, 0, 0; target = 4;
  EXPECT_EQ(target, supercell.add_offset(source, offset).first);
  source = 14; offset << 0, 1, 0; target = 15;
  EXPECT_EQ(target, supercell.add_offset(source, offset).first);
  source = 14; offset << 0, 0, 1; target = 30;
  EXPECT_EQ(target, supercell.add_offset(source, offset).first);
  source = 14; offset << -1, 0, 0; target = 13;
  EXPECT_EQ(target, supercell.add_offset(source, offset).first);
  source = 14; offset << 0, -1, 0; target = 11;
  EXPECT_EQ(target, supercell.add_offset(source, offset).first);
  source = 14; offset << 0, 0, -1; target = 62;
  EXPECT_EQ(target, supercell.add_offset(source, offset).first);
}
