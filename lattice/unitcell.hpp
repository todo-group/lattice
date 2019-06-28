/*****************************************************************************
*
* Copyright (C) 2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef LATTICE_UNITCELL_HPP
#define LATTICE_UNITCELL_HPP

#include "basis.hpp"
#include "types.hpp"

namespace lattice {

class unitcell {
public:
  struct site_t {
    site_t() {}
    site_t(const coordinate_t& pos, int tp) :
      coordinate(pos), type(tp) /* , neighbors(0), neighbor_bonds(0) */ {}
    coordinate_t coordinate;
    int type;
    // std::vector<std::size_t> neighbors;
    // std::vector<std::size_t> neighbor_bonds;
  };
  
  struct bond_t {
    bond_t() {}
    bond_t(std::size_t s, std::size_t t, offset_t os, int tp) :
      source(s), target(t), target_offset(os), type(tp) {}
    std::size_t source, target;
    offset_t target_offset;
    int type;
  };

  unitcell() {}
  explicit unitcell(std::size_t dim) : basis_(dim) {
    coordinate_t pos = coordinate_t::Zero(dim);
    add_site(pos, 0);
    for (std::size_t m = 0; m < dim; ++m) {
      offset_t os = offset_t::Zero(dim);
      os(m) = 1;
      add_bond(0, 0, os, 0);
    }
  }
  explicit unitcell(const basis& bs) : basis_(bs) {}
  unitcell(const unitcell& cell, const extent_t& extent);
  unitcell(const unitcell& cell, const span_t& span);

  std::size_t dimension() const { return basis_.dimension(); }
  double volume() const { return basis_.volume(); }
  basis_t basis_vectors() const { return basis_.basis_vectors(); }
  std::size_t num_sites() const { return sites_.size(); }
  std::size_t num_bonds() const { return bonds_.size(); }
  const site_t& site(std::size_t s) const { return sites_[s]; }
  const bond_t& bond(std::size_t b) const { return bonds_[b]; }
  std::size_t max_neighbors() const {
    std::vector<std::size_t> num_neighbors(num_sites(), 0);
    for (std::size_t b = 0; b < num_bonds(); ++b) {
      num_neighbors[bonds_[b].source] += 1;
      num_neighbors[bonds_[b].target] += 1;
    }
    return *std::max_element(num_neighbors.begin(), num_neighbors.end());
  }
  
  std::size_t add_site(const coordinate_t& pos, int tp) {
    if (pos.size() != dimension())
      throw std::invalid_argument("site coordinate dimension mismatch");
    for (std::size_t i = 0; i < dimension(); ++i) {
      if (pos[i] < 0 || pos[i] >= 1.0)
        throw std::invalid_argument("site coordinate out of range");
    }
    std::size_t s = sites_.size();
    sites_.push_back(site_t(pos, tp));
    return s;
  }
    
  std::size_t add_bond(std::size_t s, std::size_t t, const offset_t& os, int tp) {
    if (s >= num_sites() || t >= num_sites())
      throw std::invalid_argument("site index out of range");
    if (os.size() != dimension())
      throw std::invalid_argument("unitcell offset dimension mismatch");
    std::size_t b = bonds_.size();
    bonds_.push_back(bond_t(s, t, os, tp));
    // sites_[s].neighbor_bonds.push_back(b);
    // sites_[t].neighbor_bonds.push_back(b);
    // sites_[s].neighbors.push_back(t);
    // sites_[t].neighbors.push_back(s);
    return b;
  }

private:
  basis basis_;
  std::vector<site_t> sites_;
  std::vector<bond_t> bonds_;
};

std::size_t dimension(const unitcell& cell) {
  return cell.dimension();
}

} // end namespace lattice

#endif
