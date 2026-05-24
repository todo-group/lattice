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
#include "lattice/extent.hpp"

using namespace lattice;

std::size_t volume(const span_t& span) {
  return std::size_t(span.cast<double>().determinant() + 0.1);
}

TEST(ExtentTest, Extent1) {
  auto span = extent(3);
  EXPECT_EQ(3, volume(span));
  EXPECT_EQ(3, span.trace());
}

TEST(ExtentTest, Extent2) {
  auto span = extent(3, 5);
  EXPECT_EQ(15, volume(span));
  EXPECT_EQ(8, span.trace());
}

TEST(ExtentTest, Extent3) {
  auto span = extent(3, 2, 4);
  EXPECT_EQ(24, volume(span));
  EXPECT_EQ(9, span.trace());
}

TEST(ExtentTest, Extent4) {
  auto span = extent(3, 2, 4, 5);
  EXPECT_EQ(120, volume(span));
  EXPECT_EQ(14, span.trace());
}

TEST(ExtentTest, Extent5) {
  auto span = extent(3, 2, 4, 5, 1);
  EXPECT_EQ(120, volume(span));
  EXPECT_EQ(15, span.trace());
}
