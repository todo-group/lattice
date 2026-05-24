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

#ifndef LATTICE_COORDINATE_HPP
#define LATTICE_COORDINATE_HPP

#include "types.hpp"

namespace lattice {

inline coordinate_t coordinate(double x) {
  coordinate_t pos(1);
  pos << x;
  return pos;
}

inline coordinate_t coordinate(double x, double y) {
  coordinate_t pos(2);
  pos << x, y;
  return pos;
}

inline coordinate_t coordinate(double x, double y, double z) {
  coordinate_t pos(3);
  pos << x, y, z;
  return pos;
}

inline coordinate_t coordinate(double x, double y, double z, double w) {
  coordinate_t pos(4);
  pos << x, y, z, w;
  return pos;
}

inline coordinate_t coordinate(double x, double y, double z, double w, double v) {
  coordinate_t pos(5);
  pos << x, y, z, w, v;
  return pos;
}

} // end namespace lattice

#endif
