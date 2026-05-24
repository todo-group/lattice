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
#include "lattice/supercell.hpp"

using namespace lattice;

void check(const supercell& supercell) {
  for (std::size_t i = 0; i < supercell.num_cells(); ++i) {
    EXPECT_EQ(i, supercell.lcord2index(supercell.index2lcord(i)));
  }
  for (std::size_t i = 0; i < supercell.num_cells(); ++i) {
    auto offset = supercell.lcord2offset(supercell.index2lcord(i));
    //// std::cout << i << ' ' << supercell.index2lcord(i) << " (" << offset.transpose() << ") " << supercell.offset2lcord(offset) << std::endl;
    EXPECT_EQ(i, supercell.lcord2index(supercell.offset2lcord(offset)));
  }
}

TEST(SupercellTest, TwoD1) {
  span_t span(2, 2);
  span << 4, 2,
          1, 3;
  supercell supercell(span);
  offset_t offset(2);

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

  check(supercell);

  std::size_t source, target;
  offset_t crossing(2), crossing_r(2);
  source = 0; offset << 1, 0; target = 7; crossing << 0, -1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 3; offset << -1, -1; target = 8; crossing << 0, -1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 1; offset << 4, 0; target = 7; crossing << 1, -1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 7; offset << -2, 0; target = 3; crossing << -1, 1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
}

TEST(SupercellTest, TwoD2) {
  span_t span(2, 2);
  span << 4, -1,
          1, 3;
  supercell supercell(span);
  offset_t offset(2);

  EXPECT_EQ(13, supercell.num_cells());
  
  // within_supercell
  offset << 0, 0; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, 1; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 1; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, 2; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 2; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, 3; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 3; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, -1; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 0; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 1, 0; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 1; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 4, 1; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 2; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 4, 2; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 3; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 4, 3; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 0, 4; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 3, 4; EXPECT_FALSE(supercell.within_supercell(offset));

  check(supercell);

  std::size_t source, target;
  offset_t crossing(2), crossing_r(2);
  source = 0; offset << 1, 0; target = 9; crossing << 0, -1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 4; offset << -1, -1; target = 10; crossing << 0, -1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 1; offset << 4, 0; target = 0; crossing << 1, 0;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 7; offset << 3, 1; target = 6; crossing << 1, 0;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 7; offset << 1, 3; target = 1; crossing << 1, 1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
}

TEST(SupercellTest, TwoD3) {
  span_t span(2, 2);
  span << -1, 4,
           3, 1;
  supercell supercell(span);
  offset_t offset(2);

  EXPECT_EQ(13, supercell.num_cells());
  
  // within_supercell
  offset << 0, 0; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, 1; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 1; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, 2; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 2; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, 3; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 3; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, -1; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 0; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 1, 0; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 1; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 4, 1; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 2; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 4, 2; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 3; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 4, 3; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 0, 4; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 3, 4; EXPECT_FALSE(supercell.within_supercell(offset));

  check(supercell);

  std::size_t source, target;
  offset_t crossing(2), crossing_r(2);
  source = 0; offset << 1, 0; target = 9; crossing << -1, 0;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 4; offset << -1, -1; target = 10; crossing << -1, 0;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 1; offset << 4, 0; target = 0; crossing << 0, 1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 7; offset << 3, 1; target = 6; crossing << 0, 1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 7; offset << 1, 3; target = 1; crossing << 1, 1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
}

TEST(SupercellTest, TwoD4) {
  span_t span(2, 2);
  span << 4, -1,
          1,  4;
  supercell supercell(span);
  offset_t offset(2);

  EXPECT_EQ(17, supercell.num_cells());
  
  // within_supercell
  offset << 0, 0; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, 1; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 1; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, 2; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 2; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, 4; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 4; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, -1; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 0; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 1, 0; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 1; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 4, 1; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 2; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 4, 2; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 4; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 4, 4; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 0, 5; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 3, 5; EXPECT_FALSE(supercell.within_supercell(offset));

  check(supercell);

  std::size_t source, target;
  offset_t crossing(2), crossing_r(2);
  source = 0; offset << 1, 0; target = 13; crossing << 0, -1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 4; offset << -1, -1; target = 14; crossing << 0, -1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 1; offset << 4, 0; target = 0; crossing << 1, 0;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 7; offset << 3, 1; target = 6; crossing << 1, 0;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 7; offset << 1, 3; target = 0; crossing << 1, 1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
}

TEST(SupercellTest, TwoDSimple) {
  supercell supercell(extent(4, 4));
  offset_t offset(2);

  // within_supercell
  offset << 0, 0; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, 3; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 1, 0; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 1, 3; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 0; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 3, 3; EXPECT_TRUE(supercell.within_supercell(offset));
  offset << 0, 4; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 4, 4; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << 0, -1; EXPECT_FALSE(supercell.within_supercell(offset));
  offset << -1, 0; EXPECT_FALSE(supercell.within_supercell(offset));

  check(supercell);

  std::size_t source, target;
  offset_t crossing(2), crossing_r(2);
  source = 0; offset << 1, 0; target = 1; crossing << 0, 0;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 3; offset << 1, 0; target = 0; crossing << 1, 0;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 0; offset << 4, 0; target = 0; crossing << 1, 0;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
  source = 1; offset << 0, 4; target = 1; crossing << 0, 1;
  EXPECT_EQ(std::make_pair(target, crossing), supercell.add_offset(source, offset));
}

TEST(SupercellTest, ThreeD) {
  span_t span(3, 3);
  span << 4, 0, 0,
          2, 4, 0,
          0, 0, 4;
  supercell supercell(span);
  offset_t offset(3);

  check(supercell);

  std::size_t source, target;
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
