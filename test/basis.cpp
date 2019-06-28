/*****************************************************************************
*
* Copyright (C) 2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include "gtest/gtest.h"
#include "lattice/basis.hpp"

TEST(BasisTest, SimpleBasis0) {
  lattice::basis basis(0);
  std::cout << basis.dimension() << std::endl;
  std::cout << basis.basis_vectors() << std::endl;
  EXPECT_EQ(0, basis.dimension());
}

TEST(BasisTest, SimpleBasis2) {
  lattice::basis_t bs(2, 2);
  bs << 1, 0,
        0, 1.5;
  lattice::basis basis(bs);
  std::cout << basis.dimension() << std::endl;
  std::cout << basis.basis_vectors() << std::endl;
  std::cout << basis.volume() << std::endl;
  EXPECT_EQ(2, basis.dimension());
  EXPECT_DOUBLE_EQ(1.5, basis.volume());
}

TEST(BasisTest, SimpleBasis3) {
  lattice::basis basis(3);
  std::cout << basis.dimension() << std::endl;
  std::cout << basis.basis_vectors() << std::endl;
  std::cout << basis.volume() << std::endl;
  EXPECT_EQ(3, basis.dimension());
  EXPECT_DOUBLE_EQ(1, basis.volume());
}
