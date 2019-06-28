/*****************************************************************************
*
* Copyright (C) 2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include "gtest/gtest.h"
#include "lattice/lattice.hpp"

TEST(LatticeTest, SimpleSquare) {
  lattice::lattice lat(2, 4);
  std::cout << "dimension: " << lat.dimension() << std::endl
            << "number of sites: " << lat.num_sites() << std::endl
            << "number of bonds: " << lat.num_bonds() << std::endl;
  for (std::size_t s = 0; s < lat.num_sites(); ++s) {
    std::cout << "site: " << s << ' ' << lat.site_type(s) << ' '
              << "[ " << lat.coordinate(s).transpose() << " ] neighbors[ ";
    for (std::size_t k = 0; k < lat.num_neighbors(s); ++k)
      std::cout << lat.neighbor(s, k) << ' ';
    std::cout << "] neighbor_bonds[ ";
    for (std::size_t k = 0; k < lat.num_neighbors(s); ++k)
      std::cout << lat.neighbor_bond(s, k) << ' ';
    std::cout << "]\n";
  }
  
  for (std::size_t s = 0; s < lat.num_sites(); ++s) {
    for (unsigned int k = 0; k < lat.num_neighbors(k); ++k) {
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
  
  std::vector<std::size_t> count(lat.num_sites(), 0);
  for (std::size_t b = 0; b < lat.num_bonds(); ++b) {
    ++count[lat.source(b)];
    ++count[lat.target(b)];
  }
  for (std::size_t s = 0; s < lat.num_sites(); ++s) {
    EXPECT_EQ(lat.num_neighbors(s), count[s]);
  }
}

TEST(LatticeTest, SimpleCubic) {
  lattice::lattice lat(3, 3);
  std::cout << "dimension: " << lat.dimension() << std::endl
            << "number of sites: " << lat.num_sites() << std::endl
            << "number of bonds: " << lat.num_bonds() << std::endl;
  for (std::size_t s = 0; s < lat.num_sites(); ++s) {
    std::cout << "site: " << s << ' ' << lat.site_type(s) << ' '
              << "[ " << lat.coordinate(s).transpose() << " ] neighbors[ ";
    for (std::size_t k = 0; k < lat.num_neighbors(s); ++k)
      std::cout << lat.neighbor(s, k) << ' ';
    std::cout << "] neighbor_bonds[ ";
    for (std::size_t k = 0; k < lat.num_neighbors(s); ++k)
      std::cout << lat.neighbor_bond(s, k) << ' ';
    std::cout << "]\n";
  }
  
  for (std::size_t s = 0; s < lat.num_sites(); ++s) {
    for (unsigned int k = 0; k < lat.num_neighbors(k); ++k) {
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
  
  std::vector<std::size_t> count(lat.num_sites(), 0);
  for (std::size_t b = 0; b < lat.num_bonds(); ++b) {
    ++count[lat.source(b)];
    ++count[lat.target(b)];
  }
  for (std::size_t s = 0; s < lat.num_sites(); ++s) {
    EXPECT_EQ(lat.num_neighbors(s), count[s]);
  }
}
