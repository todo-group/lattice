/*****************************************************************************
*
* Copyright (C) 2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef LATTICE_LATTICE_HPP
#define LATTICE_LATTICE_HPP

#include "supercell.hpp"
#include "unitcell.hpp"

namespace lattice {

class lattice {
public:
  lattice() {}
  lattice(std::size_t dim, std::size_t length) :
    unitcell_(dim), supercell_(dim, length), boundary_(dim, boundary_t::periodic) {
    init();
  }
  lattice(std::size_t dim, const extent_t& extent) :
    unitcell_(dim), supercell_(extent), boundary_(dim, boundary_t::periodic) {
    init();
  }
  lattice(const unitcell& cell, std::size_t length) :
    unitcell_(cell), supercell_(cell.dimension(), length),
    boundary_(cell.dimension(), boundary_t::periodic) {
    init();
  }
  lattice(const unitcell& cell, const extent_t& extent) :
    unitcell_(cell), supercell_(extent), boundary_(cell.dimension(), boundary_t::periodic) {
    init();
  }
  lattice(const unitcell& cell, const extent_t& extent, const std::vector<boundary_t>& boundary) :
    unitcell_(cell), supercell_(extent), boundary_(boundary) {
    init();
  }
  lattice(const unitcell& cell, const span_t& span, const std::vector<boundary_t>& boundary) :
    unitcell_(cell), supercell_(span), boundary_(boundary) {
    init();
  }

  void init() {
    if (unitcell_.dimension() != supercell_.dimension() ||
        unitcell_.dimension() != boundary_.size())
      throw std::invalid_argument("dimension mismatch");
    std::size_t dim = unitcell_.dimension();
    std::size_t ns = unitcell_.num_sites() * supercell_.num_cells();
    std::size_t nb = unitcell_.num_bonds() * supercell_.num_cells();
    std::size_t max_neighbors = unitcell_.max_neighbors();
    
    site_types_.resize(ns);
    site_coordinates_.resize(dim, ns);
    num_neighbors_.resize(ns); std::fill(num_neighbors_.begin(), num_neighbors_.end(), 0);
    neighbors_.resize(ns, max_neighbors);
    neighbor_bonds_.resize(ns, max_neighbors);
    bond_types_.resize(nb);
    edge_sites_.resize(nb);
    
    std::size_t s = 0;
    for (std::size_t c = 0; c < supercell_.num_cells(); ++c) {
      auto cell_offset = supercell_.offset(c);
      for (std::size_t t = 0; t < unitcell_.num_sites(); ++t) {
        site_types_[s] = unitcell_.site(t).type;
        site_coordinates_.col(s) = unitcell_.basis_vectors() *
          (cell_offset.cast<double>() + unitcell_.site(t).coordinate);
        ++s;
      }
    }
    
    std::size_t b = 0;
    for (std::size_t c = 0; c < supercell_.num_cells(); ++c) {
      for (std::size_t u = 0; u < unitcell_.num_bonds(); ++u) {
        bond_types_[b] = unitcell_.bond(u).type;
        std::size_t s = c * unitcell_.num_sites() + unitcell_.bond(u).source;
        std::size_t t = supercell_.add_offset(c, unitcell_.bond(u).target_offset).first *
          unitcell_.num_sites() + unitcell_.bond(u).target;
        neighbors_(s, num_neighbors_[s]) = t;
        neighbors_(t, num_neighbors_[t]) = s;
        neighbor_bonds_(s, num_neighbors_[s]) = b;
        neighbor_bonds_(t, num_neighbors_[t]) = b;
        ++num_neighbors_[s];
        ++num_neighbors_[t];
        edge_sites_[b] = std::make_pair(s, t);
        ++b;
      }
    }
  }
  
  std::size_t dimension() const { return unitcell_.dimension(); }
  
  std::size_t num_sites() const { return site_types_.size(); }
  std::size_t site_type(std::size_t s) const { return site_types_[s]; }
  coordinate_t coordinate(std::size_t s) const { return site_coordinates_.col(s); }
  double coordinate(std::size_t s, std::size_t d) const { return site_coordinates_(s, d); }
  std::size_t num_neighbors(std::size_t s) const { return num_neighbors_[s]; }
  std::size_t neighbor(std::size_t s, std::size_t k) const { return neighbors_(s, k); }
  std::size_t neighbor_bond(std::size_t s, std::size_t k) const { return neighbor_bonds_(s, k); }

  std::size_t num_bonds() const { return bond_types_.size(); }
  std::size_t bond_type(std::size_t b) const { return bond_types_[b]; }
  std::size_t source(std::size_t b) const { return edge_sites_[b].first; }
  std::size_t target(std::size_t b) const { return edge_sites_[b].second; }
  const std::pair<std::size_t, std::size_t>& edge_sites(std::size_t b) const {
    return edge_sites_[b];
  }
      
private:
  unitcell unitcell_;
  supercell supercell_;
  std::vector<boundary_t> boundary_;
  
  std::vector<int> site_types_;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> site_coordinates_;
  std::vector<std::size_t> num_neighbors_;
  Eigen::Matrix<std::size_t, Eigen::Dynamic, Eigen::Dynamic> neighbors_;
  Eigen::Matrix<std::size_t, Eigen::Dynamic, Eigen::Dynamic> neighbor_bonds_;

  std::vector<int> bond_types_;
  std::vector<std::pair<std::size_t, std::size_t>> edge_sites_;
};

} // end namespace lattice

#endif // ALPS_LATTICE_LATTICE_HPP
