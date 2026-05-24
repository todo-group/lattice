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
#include "lattice/basis.hpp"

using namespace lattice;

TEST(BasisTest, SimpleBasis0) {
  basis basis = basis::simple(0);
  std::cout << basis.dimension() << std::endl;
  std::cout << basis.basis_vectors() << std::endl;
  EXPECT_EQ(0, basis.dimension());
}

TEST(BasisTest, SimpleBasis2) {
  basis_t bs(2, 2);
  bs << 1, 0,
        0, 1.5;
  basis basis(bs);
  std::cout << basis.dimension() << std::endl;
  std::cout << basis.basis_vectors() << std::endl;
  std::cout << basis.volume() << std::endl;
  EXPECT_EQ(2, basis.dimension());
  EXPECT_DOUBLE_EQ(1.5, basis.volume());
}

TEST(BasisTest, SimpleBasis3) {
  basis basis = basis::simple(3);
  std::cout << basis.dimension() << std::endl;
  std::cout << basis.basis_vectors() << std::endl;
  std::cout << basis.volume() << std::endl;
  EXPECT_EQ(3, basis.dimension());
  EXPECT_DOUBLE_EQ(1, basis.volume());
}
