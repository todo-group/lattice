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

#ifndef LATTICE_LATTICE_HPP
#define LATTICE_LATTICE_HPP

#include <vector>
#include "supercell.hpp"
#include "unitcell.hpp"

namespace lattice {

class lattice {
public:
  struct site_t {
    site_t() {}
    site_t(int tp) : type(tp), neighbors(0), neighbor_bonds(0) {}
    int type;
    std::vector<std::size_t> neighbors;
    std::vector<std::size_t> neighbor_bonds;
  };
  
  struct bond_t {
    bond_t() {}
    bond_t(std::size_t s, std::size_t t, int tp) : source(s), target(t), type(tp) {}
    std::size_t source, target;
    int type;
  };

  lattice() {}
  lattice(std::size_t dim, std::size_t length) {
    init(unitcell(dim), supercell(dim, length), std::vector<boundary_t>(dim, boundary_t::periodic));
  }
  lattice(std::size_t dim, const extent_t& extent) {
    init(unitcell(dim), supercell(extent), std::vector<boundary_t>(dim, boundary_t::periodic));
  }
  lattice(const unitcell& cell, std::size_t length) {
    init(cell, supercell(cell.dimension(), length),
         std::vector<boundary_t>(cell.dimension(), boundary_t::periodic));
  }
  lattice(const unitcell& cell, const extent_t& extent) {
    init(cell, supercell(extent), std::vector<boundary_t>(cell.dimension(), boundary_t::periodic));
  }
  lattice(const unitcell& cell, const extent_t& extent, const std::vector<boundary_t>& boundary) {
    init(cell, supercell(extent), boundary);
  }
  lattice(const unitcell& cell, const span_t& span, const std::vector<boundary_t>& boundary) {
    init(cell, supercell(span), boundary);
  }

  void init(const unitcell& cell, const supercell& super, const std::vector<boundary_t>& boundary) {
    if (cell.dimension() != super.dimension() || cell.dimension() != boundary.size())
      throw std::invalid_argument("dimension mismatch");

    dim_ = cell.dimension();
    sites_.clear();
    coordinates_.clear();
    bonds_.clear();

    for (std::size_t c = 0; c < super.num_cells(); ++c) {
      auto cell_offset = super.offset(c);
      for (std::size_t t = 0; t < cell.num_sites(); ++t) {
        coordinate_t pos = cell.basis_vectors() *
          (cell_offset.cast<double>() + cell.site(t).coordinate);
        add_site(cell.site(t).type, pos);
      }
    }
    for (std::size_t c = 0; c < super.num_cells(); ++c) {
      for (std::size_t u = 0; u < cell.num_bonds(); ++u) {
        std::size_t b = bonds_.size();
        std::size_t s = c * cell.num_sites() + cell.bond(u).source;
        std::size_t t = super.add_offset(c, cell.bond(u).target_offset).first *
          cell.num_sites() + cell.bond(u).target;
        add_bond(s, t, cell.bond(u).type);
      }
    }
  }

  std::size_t add_site(int tp, const coordinate_t& pos) {
    if (pos.size() != dim_)
      throw std::invalid_argument("dimension mismatch");
    std::size_t s = sites_.size();
    sites_.push_back(tp);
    coordinates_.push_back(pos);
    return s;
  }

  std::size_t add_bond(std::size_t s, std::size_t t, int tp) {
    if (s >= sites_.size() || t >= sites_.size())
      throw std::invalid_argument("site index out of range");
    std::size_t b = bonds_.size();
    bonds_.push_back(bond_t(s, t, tp));
    sites_[s].neighbors.push_back(t);
    sites_[s].neighbor_bonds.push_back(b);
    sites_[t].neighbors.push_back(s);
    sites_[t].neighbor_bonds.push_back(b);
    return b;
  }
      
  std::size_t dimension() const { return dim_; }
  
  std::size_t num_sites() const { return sites_.size(); }
  std::size_t site_type(std::size_t s) const { return sites_[s].type; }
  const coordinate_t& coordinate(std::size_t s) const { return coordinates_[s]; }
  std::size_t num_neighbors(std::size_t s) const {
    return sites_[s].neighbors.size();
  }
  std::size_t neighbor(std::size_t s, std::size_t k) const {
    return sites_[s].neighbors[k];
  }
  std::size_t neighbor_bond(std::size_t s, std::size_t k) const {
    return sites_[s].neighbor_bonds[k];
  }

  std::size_t num_bonds() const { return bonds_.size(); }
  std::size_t bond_type(std::size_t b) const { return bonds_[b].type; }
  std::size_t source(std::size_t b) const { return bonds_[b].source; }
  std::size_t target(std::size_t b) const { return bonds_[b].target; }
  std::pair<std::size_t, std::size_t> edge_sites(std::size_t b) const {
    return std::make_pair(bonds_[b].source, bonds_[b].target);
  }
      
private:
  std::size_t dim_;
  std::vector<site_t> sites_;
  std::vector<coordinate_t> coordinates_;
  std::vector<bond_t> bonds_;
};

} // end namespace lattice

#endif // LATTICE_LATTICE_HPP
