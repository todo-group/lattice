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

#ifndef LATTICE_OFFSET_HPP
#define LATTICE_OFFSET_HPP

#include "types.hpp"

namespace lattice {

inline offset_t offset(long x) {
  offset_t os(1);
  os << x;
  return os;
}

inline offset_t offset(long x, long y) {
  offset_t os(2);
  os << x, y;
  return os;
}

inline offset_t offset(long x, long y, long z) {
  offset_t os(3);
  os << x, y, z;
  return os;
}

inline offset_t offset(long x, long y, long z, long w) {
  offset_t os(4);
  os << x, y, z, w;
  return os;
}

inline offset_t offset(long x, long y, long z, long w, long v) {
  offset_t os(5);
  os << x, y, z, w, v;
  return os;
}

} // end namespace lattice

#endif
