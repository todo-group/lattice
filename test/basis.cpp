#include "gtest/gtest.h"
#include "lattice/basis.hpp"

TEST(BasisTest, SimpleBasis0) {
  std::string name = "no dimension";
  lattice::basis basis(name, 0);
  std::cout << basis.name() << std::endl;
  std::cout << basis.dimension() << std::endl;
  std::cout << basis.basis_vectors() << std::endl;
  EXPECT_EQ(0, basis.dimension());
}

TEST(BasisTest, SimpleBasis2) {
  std::string name = "square lattice";
  lattice::basis_t bs(2, 2);
  bs << 1, 0,
        0, 1.5;
  lattice::basis basis(name, bs);
  std::cout << basis.name() << std::endl;
  std::cout << basis.dimension() << std::endl;
  std::cout << basis.basis_vectors() << std::endl;
  std::cout << basis.volume() << std::endl;
  EXPECT_EQ(2, basis.dimension());
  EXPECT_DOUBLE_EQ(1.5, basis.volume());
}

TEST(BasisTest, SimpleBasis3) {
  lattice::basis basis("simple cubic lattice", 3);
  std::cout << basis.name() << std::endl;
  std::cout << basis.dimension() << std::endl;
  std::cout << basis.basis_vectors() << std::endl;
  std::cout << basis.volume() << std::endl;
  EXPECT_EQ(3, basis.dimension());
  EXPECT_DOUBLE_EQ(1, basis.volume());
}
