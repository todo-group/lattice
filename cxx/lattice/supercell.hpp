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

#ifndef LATTICE_SUPERCELL_HPP
#define LATTICE_SUPERCELL_HPP

#include <algorithm>
#include <exception>
#include <vector>
#include "extent.hpp"
#include "types.hpp"

namespace lattice {

class supercell {
public:
  supercell() {}
  supercell(std::size_t dim, std::size_t length) {
    span_t span = length * span_t::Identity(dim, dim);
    init(span);
  }
  supercell(const span_t& span) { init(span); }

  void init(const span_t& span) {
    span_ = span;
    if (span_.cols() != span_.rows())
      throw std::invalid_argument("invalid shape of span matrix");
    dim_ = span_.cols();
    num_cells_ = std::size_t(std::abs(span_.cast<double>().determinant()) + 0.5);
    rs_ = span_.cast<double>().inverse();

    nmin_ = offset_t::Zero(dim_);
    nmax_ = offset_t::Zero(dim_);
    for (std::size_t i = 0; i < std::size_t(1 << dim_); ++i) {
      for (std::size_t m = 0; m < dim_; ++m) {
        long v = 0;
        for (std::size_t n = 0; n < dim_; ++n) v += ((i >> n) & 1) * span_(m, n);
        nmin_(m) = std::min(nmin_(m), v);
        nmax_(m) = std::max(nmax_(m), v);
      }
    }
    
    std::size_t lcord_max = 1;
    for (std::size_t m = 0; m < dim_; ++m) lcord_max *= (nmax_(m) - nmin_(m));
    lcord2index_.resize(lcord_max);
    index2lcord_.clear();
    for (std::size_t i = 0; i < lcord2index_.size(); ++i) {
      offset_t c = lcord2offset(i);
      if (within_supercell(c)) {
        lcord2index_[i] = index2lcord_.size();
        index2lcord_.push_back(i);
      }
    }
    if (num_cells_ != index2lcord_.size())
      throw std::logic_error("supercell::init() internal error");
  }

  std::size_t dimension() const { return dim_; }
  std::size_t num_cells() const { return num_cells_; }

  std::pair<std::size_t, offset_t> add_offset(std::size_t index, const offset_t& offset) const {
    offset_t cell = lcord2offset(index2lcord_.at(index)) + offset;
    offset_t crossing = offset_t::Zero(dim_);
    const double eps = 1.0e-8;
    bool checked = false;
    while (!checked) {
      checked = true;
      auto p = rs_ * cell.cast<double>();
      for (std::size_t m = 0; m < dim_; ++m) {
        if (p(m) < -eps) {
          checked = false;
          cell += span_.col(m);
          crossing(m) -= 1;
        } else if (p(m) > (1.0 - eps)) {
          checked = false;
          cell -= span_.col(m);
          crossing(m) += 1;
        }
      }
    }
    return std::make_pair(lcord2index(offset2lcord(cell)), crossing);
  }

  offset_t offset(std::size_t index) const {
    return lcord2offset(index2lcord(index));
  }

  std::size_t lcord2index(std::size_t lc) const { return lcord2index_.at(lc); }

  std::size_t index2lcord(std::size_t index) const { return index2lcord_.at(index); }

  offset_t lcord2offset(std::size_t lc) const {
    offset_t c(dim_);
    for (std::size_t m = 0; m < dim_; ++m) {
      c(m) = (lc % (nmax_(m) - nmin_(m))) + nmin_(m);
      lc /= (nmax_(m) - nmin_(m));
    }
    return c;
  }
  
  std::size_t offset2lcord(const offset_t& offset) const {
    std::size_t lc = 0;
    for (std::size_t m = dim_; m > 0; --m) {
      lc *= (nmax_(m-1) - nmin_(m-1));
      lc += (offset(m-1) - nmin_(m-1));
    }
    return lc;
  }

  bool within_supercell(const offset_t& cell) const {
    const double eps = 1.0e-8;
    auto p = rs_ * cell.cast<double>();
    bool res = true;
    for (std::size_t m = 0; m < dim_; ++m) {
      if (p(m) < -eps || p(m) > (1.0 - eps)) res = false;
    }
    return res;
  }
 
private:
  std::size_t dim_, num_cells_;
  span_t span_;
  Eigen::MatrixXd rs_;
  offset_t nmin_, nmax_;
  std::vector<std::size_t> lcord2index_, index2lcord_;
};

} // end namespace lattice

#endif // LATTICE_SUPERCELL_HPP
