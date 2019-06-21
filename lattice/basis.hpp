#ifndef LATTICE_BASIS_HPP
#define LATTICE_BASIS_HPP

#include <exception>
#include "types.hpp"

namespace lattice {

class basis {
public:
  basis() {}
  basis(const std::string& name, const basis_t& bs) : name_(name) { set_basis(bs); }
  basis(const std::string& name, std::size_t dim) : name_(name) {
    set_basis(basis_t::Identity(dim, dim));
  }
  void clear() { *this = basis(); }
  
  void set_name(const std::string& name) { name_ = name; }
  void set_basis(const basis_t& bs) {
    basis_ = bs;
    if (basis_.rows() != basis_.cols())
      throw std::runtime_error("basis dimension mismatch");
  }
  const std::string& name() const { return name_; }
  std::size_t dimension() const { return basis_.rows(); }
  basis_t basis_vectors() const { return basis_; }
private:
  std::string name_;
  basis_t basis_;
};

std::size_t dimension(const basis& bs) {
  return bs.dimension();
}
  
} // end namespace lattice

#endif
