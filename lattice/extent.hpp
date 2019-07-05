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

#ifndef LATTICE_EXTENT_HPP
#define LATTICE_EXTENT_HPP

#include "types.hpp"

namespace lattice {

inline span_t extent(std::size_t L1) {
  span_t span(1, 1); span.setZero();
  span(0, 0) = L1;
  return span;
}

inline span_t extent(std::size_t L1, std::size_t L2) {
  span_t span(2, 2); span.setZero();
  span(0, 0) = L1;
  span(1, 1) = L2;
  return span;
}

inline span_t extent(std::size_t L1, std::size_t L2, std::size_t L3) {
  span_t span(3, 3); span.setZero();
  span(0, 0) = L1;
  span(1, 1) = L2;
  span(2, 2) = L3;
  return span;
}

inline span_t extent(std::size_t L1, std::size_t L2, std::size_t L3, std::size_t L4) {
  span_t span(4, 4); span.setZero();
  span(0, 0) = L1;
  span(1, 1) = L2;
  span(2, 2) = L3;
  span(3, 3) = L4;
  return span;
}

inline span_t extent(std::size_t L1, std::size_t L2, std::size_t L3, std::size_t L4,
                     std::size_t L5) {
  span_t span(5, 5); span.setZero();
  span(0, 0) = L1;
  span(1, 1) = L2;
  span(2, 2) = L3;
  span(3, 3) = L4;
  span(4, 4) = L5;
  return span;
}

} // end namespace lattice

#endif
