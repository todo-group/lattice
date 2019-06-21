#include "gtest/gtest.h"
#include "lattice/basis.hpp"

TEST(BasisTest, SimpleBasis2) {
  std::string name = "square lattice";
  lattice::basis_t bs(2, 2);
  bs << 1, 0,
        0, 1;
  lattice::basis basis(name, bs);
  std::cout << basis.name() << std::endl;
  std::cout << basis.dimension() << std::endl;
  std::cout << basis.basis_vectors() << std::endl;
}

TEST(BasisTest, SimpleBasis3) {
  lattice::basis basis("simple cubic lattice", 3);
  std::cout << basis.name() << std::endl;
  std::cout << basis.dimension() << std::endl;
  std::cout << basis.basis_vectors() << std::endl;
}
