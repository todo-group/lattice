/*****************************************************************************
*
* Copyright (C) 2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef LATTICE_BASIS_HPP
#define LATTICE_BASIS_HPP

#include <exception>
#include "types.hpp"

namespace lattice {

class basis {
public:
  basis() {}
  explicit basis(const basis_t& bs) { set_basis(bs); }
  explicit basis(std::size_t dim) { set_basis(basis_t::Identity(dim, dim)); }
  void clear() { *this = basis(); }
  
  void set_basis(const basis_t& bs) {
    basis_ = bs;
    if (basis_.rows() != basis_.cols())
      throw std::runtime_error("basis dimension mismatch");
  }
  std::size_t dimension() const { return basis_.rows(); }
  basis_t basis_vectors() const { return basis_; }
  double volume() const { return std::abs(basis_.determinant()); }
      
private:
  basis_t basis_;
};

std::size_t dimension(const basis& bs) {
  return bs.dimension();
}
  
} // end namespace lattice

#endif
