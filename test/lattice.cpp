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
#include "lattice/lattice.hpp"

void check(const lattice::lattice& lat) {
  for (std::size_t s = 0; s < lat.num_sites(); ++s) {
    for (unsigned int k = 0; k < lat.num_neighbors(s); ++k) {
      EXPECT_TRUE(lat.source(lat.neighbor_bond(s, k)) == s ||
                  lat.target(lat.neighbor_bond(s, k)) == s);
      EXPECT_TRUE(lat.edge_sites(lat.neighbor_bond(s, k)).first == s ||
                  lat.edge_sites(lat.neighbor_bond(s, k)).second == s);
    }
  }

  for (std::size_t b = 0; b < lat.num_bonds(); ++b) {
    std::size_t s = lat.source(b);
    std::size_t t = lat.target(b);
    EXPECT_FALSE(s == t);
    bool has = false;
    for (std::size_t k = 0; k < lat.num_neighbors(s); ++k)
      has |= lat.neighbor_bond(s, k) == b;
    EXPECT_TRUE(has);
    has = false;
    for (std::size_t k = 0; k < lat.num_neighbors(t); ++k)
      has |= lat.neighbor_bond(t, k) == b;
    EXPECT_TRUE(has);
  }

  // detection of parallel edges
  Eigen::MatrixXi adj = Eigen::MatrixXi::Zero(lat.num_sites(), lat.num_sites());
  for (std::size_t b = 0; b < lat.num_bonds(); ++b) {
    EXPECT_EQ(0, adj(lat.source(b), lat.target(b)));
    EXPECT_EQ(0, adj(lat.target(b), lat.source(b)));
    adj(lat.source(b), lat.target(b)) = 1;
    adj(lat.target(b), lat.source(b)) = 1;
  }
  
  std::vector<std::size_t> count(lat.num_sites(), 0);
  for (std::size_t b = 0; b < lat.num_bonds(); ++b) {
    ++count[lat.source(b)];
    ++count[lat.target(b)];
  }
  for (std::size_t s = 0; s < lat.num_sites(); ++s) {
    EXPECT_EQ(lat.num_neighbors(s), count[s]);
  }
}

TEST(LatticeTest, Chain) {
  lattice::lattice lat = lattice::lattice::simple(1, 10);
  lat.print(std::cout);
  check(lat);
}

TEST(LatticeTest, SimpleSquare) {
  lattice::lattice lat = lattice::lattice::simple(2, 4);
  lat.print(std::cout);
  check(lat);
}

TEST(LatticeTest, SimpleCubic) {
  lattice::lattice lat = lattice::lattice::simple(3, 4);
  lat.print(std::cout);
  check(lat);
}

TEST(LatticeTest, SkewedSquare) {
  lattice::basis bs = lattice::basis::simple(2);
  lattice::unitcell cell = lattice::unitcell::simple(2);
  lattice::span_t span(2, 2);
  span << 4, -1,
          1,  3;
  lattice::lattice lat("skewed square", bs, cell, span);
  lat.print(std::cout);
  check(lat);
}

TEST(LatticeTest, SkewedSquareOpen) {
  lattice::basis bs = lattice::basis::simple(2);
  lattice::unitcell cell = lattice::unitcell::simple(2);
  lattice::span_t span(2, 2);
  span << 4, -1,
          1,  4;
  lattice::lattice lat("skewed square", bs, cell, span, lattice::boundary_t::open);
  lat.print(std::cout);
  check(lat);
}
