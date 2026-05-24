/*
   Copyright (C) 2019-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

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
#include "lattice/coloring.hpp"

using namespace lattice;

void check(const graph& lat) {
  auto color = coloring(lat);
  EXPECT_EQ(color.size(), lat.num_sites());
  for (std::size_t s = 0; s < lat.num_sites(); ++s) {
    EXPECT_TRUE(color[s] == 0 || color[s] == 1);
  }
  for (std::size_t b = 0; b < lat.num_bonds(); ++b) {
    EXPECT_TRUE(color[lat.source(b)] != color[lat.target(b)]);
  }
}

TEST(ColoringTest, Chain) {
  check(graph::simple(1, 10));
}

TEST(ColoringTest, ChainFail) {
  auto color = coloring(graph::simple(1, 9));
  EXPECT_EQ(0, color.size());
}

TEST(ColoringTest, SimpleSquare) {
  check(graph::simple(2, 4));
}

TEST(ColoringTest, SimpleCubic) {
  check(graph::simple(3, 4));
}
