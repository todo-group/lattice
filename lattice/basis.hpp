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

#ifndef LATTICE_BASIS_HPP
#define LATTICE_BASIS_HPP

#include <exception>
#include "types.hpp"

namespace lattice {

class basis {
public:
  basis() {}
  explicit basis(const basis_t& bs) : name_("") { set_basis(bs); }
  explicit basis(const basis_t& bs, const std::string& name) : name_(name) { set_basis(bs); }
  void clear() { *this = basis(); }
  
  void set_basis(const basis_t& bs) {
    basis_ = bs;
    if (basis_.rows() != basis_.cols())
      throw std::runtime_error("basis dimension mismatch");
  }
  std::size_t dimension() const { return basis_.rows(); }
  basis_t basis_vectors() const { return basis_; }
  double volume() const { return std::abs(basis_.determinant()); }
  std::string name() const { return name_; }

  static basis simple(std::size_t dim) {
    return basis(basis_t::Identity(dim, dim), "simple");
  }
  static basis simple(std::size_t dim, const std::string& name) {
    return basis(basis_t::Identity(dim, dim), name);
  }
  
private:
  basis_t basis_;
  std::string name_;
};

std::size_t dimension(const basis& bs) {
  return bs.dimension();
}
  
} // end namespace lattice

#endif
