/*
   Copyright (C) 2015-2020 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

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

#pragma once

#include <cmath>
#include <limits>

namespace standards {

// Solve f(x) = 0
//   Return: pair of solution and number of iteractions
//           iteration = 0 means that algorithm does not converge during max_iter.

template<typename F, typename T>
std::pair<T, std::size_t> newton_1d(F const& func, T x, std::size_t max_iter = 256) {
  using std::abs;
  std::size_t iter = 0;
  auto v = func(x);
  for (; iter < max_iter; ++iter) {
    auto xn = x - v.derivative(0) / v.derivative(1);
    auto vn = func(xn);
    if (abs(xn - x) < 2 * abs(x) * std::numeric_limits<T>::epsilon() ||
        abs(vn.derivative(0)) < 2 * std::numeric_limits<T>::epsilon())
      return std::make_pair(xn, iter + 1);
    x = xn;
    v = vn;
  }
  return std::make_pair(x, 0);
}

}
