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

#ifndef LATTICE_GRAPH_HPP
#define LATTICE_GRAPH_HPP

#include <iostream>
#include <exception>
#include <vector>
#include "basis.hpp"
#include "unitcell.hpp"
#include "supercell.hpp"

namespace lattice {

class graph {
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
  
  graph() : dim_(0) {}
  explicit graph(std::size_t dim) : dim_(dim) {}
  graph(const basis& bs, const unitcell& cell, std::size_t length,
        boundary_t boundary = boundary_t::periodic) {
    init(bs, cell, supercell(cell.dimension(), length),
         std::vector<boundary_t>(cell.dimension(), boundary));
  }
  graph(const basis& bs, const unitcell& cell, const extent_t& extent,
        boundary_t boundary = boundary_t::periodic) {
    init(bs, cell, supercell(extent), std::vector<boundary_t>(cell.dimension(), boundary));
  }
  graph(const basis& bs, const unitcell& cell, const extent_t& extent,
        const std::vector<boundary_t>& boundary) {
    init(bs, cell, supercell(extent), boundary);
  }
  graph(const basis& bs, const unitcell& cell, const span_t& span,
        boundary_t boundary = boundary_t::periodic) {
    init(bs, cell, supercell(span), std::vector<boundary_t>(cell.dimension(), boundary));
  }
  graph(const basis& bs, const unitcell& cell, const span_t& span,
        const std::vector<boundary_t>& boundary) {
    init(bs, cell, supercell(span), boundary);
  }
  
  void init(const basis& bs, const unitcell& cell, const supercell& super,
            const std::vector<boundary_t>& boundary) {
    if (cell.dimension() != super.dimension() || cell.dimension() != boundary.size())
      throw std::invalid_argument("dimension mismatch");

    dim_ = cell.dimension();
    sites_.clear();
    coordinates_.clear();
    bonds_.clear();

    for (std::size_t c = 0; c < super.num_cells(); ++c) {
      auto cell_offset = super.offset(c);
      for (std::size_t t = 0; t < cell.num_sites(); ++t) {
        coordinate_t pos = bs.basis_vectors() *
          (cell_offset.cast<double>() + cell.site(t).coordinate);
        add_site(pos, cell.site(t).type);
      }
    }
    for (std::size_t c = 0; c < super.num_cells(); ++c) {
      for (std::size_t u = 0; u < cell.num_bonds(); ++u) {
        std::size_t s = c * cell.num_sites() + cell.bond(u).source;
        std::size_t target_cell;
        offset_t cross;
        std::tie(target_cell, cross) = super.add_offset(c, cell.bond(u).target_offset);
        bool valid = true;
        for (std::size_t m = 0; m < dim_; ++m)
          if (boundary[m] == boundary_t::open && cross(m) != 0) valid = false;
        if (valid) {
          std::size_t t = target_cell * cell.num_sites() + cell.bond(u).target;
          if (s != t) add_bond(s, t, cell.bond(u).type);
        }
      }
    }
  }

  std::size_t add_site(const coordinate_t& pos, int tp) {
    if (std::size_t(pos.size()) != dim_)
      throw std::invalid_argument("dimension mismatch");
    std::size_t s = sites_.size();
    sites_.push_back(tp);
    coordinates_.push_back(pos);
    return s;
  }

  std::size_t add_bond(std::size_t s, std::size_t t, int tp) {
    if (s >= sites_.size() || t >= sites_.size())
      throw std::invalid_argument("site index out of range");
    if (s == t)
      throw std::invalid_argument("self loop is not allowed");
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
  int site_type(std::size_t s) const { return sites_[s].type; }
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
  int bond_type(std::size_t b) const { return bonds_[b].type; }
  std::size_t source(std::size_t b) const { return bonds_[b].source; }
  std::size_t target(std::size_t b) const { return bonds_[b].target; }
  std::pair<std::size_t, std::size_t> edge_sites(std::size_t b) const {
    return std::make_pair(bonds_[b].source, bonds_[b].target);
  }

  static graph simple(std::size_t dim, std::size_t length) {
    return graph(basis::simple(dim), unitcell::simple(dim), length);
  }

  static graph fully_connected(std::size_t num_sites) {
    graph g(0);
    coordinate_t pos(0);
    for (unsigned int s = 0; s < num_sites; ++s) g.add_site(pos, 0);
    for (unsigned int s = 0; s < num_sites; ++s) {
      for (unsigned int t = s + 1; t < num_sites; ++t) g.add_bond(s, t, 0);
    }
    return g;
  }
  
  void print(std::ostream& os = std::cout) const {
    os << "dimension: " << dimension() << std::endl
       << "number of sites: " << num_sites() << std::endl
       << "number of bonds: " << num_bonds() << std::endl;
    for (std::size_t s = 0; s < num_sites(); ++s) {
      os << "site: " << s << " type: " << site_type(s) << ' '
         << "( " << coordinate(s).transpose() << " ) neighbors[ ";
      for (std::size_t k = 0; k < num_neighbors(s); ++k)
        os << neighbor(s, k) << ' ';
      os << "] neighbor_bonds[ ";
      for (std::size_t k = 0; k < num_neighbors(s); ++k)
        os << neighbor_bond(s, k) << ' ';
      os << "]\n";
    }
  }
  
private:
  std::size_t dim_;
  std::vector<site_t> sites_;
  std::vector<coordinate_t> coordinates_;
  std::vector<bond_t> bonds_;
};

} // end namespace lattice

#endif // LATTICE_GRAPH_HPP
