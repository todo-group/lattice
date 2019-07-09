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

#ifndef LATTICE_TYPES_HPP
#define LATTICE_TYPES_HPP

#include <string>
#include <Eigen/Dense>

namespace lattice {

enum boundary_t { open, periodic };

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> basis_t;

// set of spanning vectors
typedef Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic> span_t;

// extent vector = diagonal spanning matrix
typedef Eigen::Matrix<long, Eigen::Dynamic, 1> extent_t;

// coordinate vector
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> coordinate_t;

// unitcell offset
typedef Eigen::Matrix<long, Eigen::Dynamic, 1> offset_t;

} // end namespace lattice

#endif
