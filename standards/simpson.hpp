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
#include <stdexcept>

namespace standards {

template<typename F, typename T, typename I>
auto simpson_1d(F const& func, T x0, T x1, I n) -> decltype(func(x0)) {
  typedef decltype(func(x0)) result_t;
  if (n == 0 || (n % 2) != 0)
    throw(std::invalid_argument("n should be positive and a multiple of two"));
  auto dx = (x1 - x0) / n;
  result_t g(0);
  // i == 0 || i == n
  g += func(x0) + func(x1);
  // i = 1/2 ... n-1/2
  for (I i = 0; i < n; ++i) {
    auto x = x0 + dx * (2*i + 1) / 2;
    g += 4 * func(x);
  }
  // i = 1 ... n-1
  for (I i = 1; i < n; ++i) {
    auto x = x0 + dx * i;
    g += 2 * func(x);
  }
  g *= dx / 6;
  return g;
}

template<typename F, typename T, typename I>
auto simpson_2d(F const& func, T x0, T y0, T x1, T y1, I nx, I ny) -> decltype(func(x0, y0)) {
  typedef decltype(func(x0, y0)) result_t;
  if (nx == 0 || (nx % 2) != 0 || ny == 0 || (ny % 2) != 0)
    throw(std::invalid_argument("nx and ny should be positive and multiples of two"));
  auto dx = (x1 - x0) / nx;
  auto dy = (y1 - y0) / ny;
  result_t g(0);
  // i == 0
  {
    auto x = x0;
    //   j == 0 || j == ny
    g += func(x, y0) + func(x, y1);
    //   j = 1/2 ... ny-1/2
    for (I j = 0; j < ny; ++j) {
      auto y = y0 + dy * (2*j + 1) / 2;
      g += 4 * func(x, y);
    }
    //   j = 1 ... ny-1
    for (I j = 1; j < ny; ++j) {
      auto y = y0 + dy * j;
      g += 2 * func(x, y);
    }
  }
  // i = 1/2 ... nx-1/2
  for (I i = 0; i < nx; ++i) {
    auto x = x0 + dx * (2*i + 1) / 2;
    //   j == 0 || j == ny
    g += 4 * (func(x, y0) + func(x, y1));
    //   j = 1/2 ... ny-1/2
    for (I j = 0; j < ny; ++j) {
      auto y = y0 + dy * (2*j + 1) / 2;
      g += 16 * func(x, y);
    }
    //   j = 1 ... ny-1
    for (I j = 1; j < ny; ++j) {
      auto y = y0 + dy * j;
      g += 8 * func(x, y);
    }
  }
  // i = 1 ... nx-1
  for (I i = 1; i < nx; ++i) {
    auto x = x0 + dx * i;
    //   j == 0 || j == ny
    g += 2 * (func(x, y0) + func(x, y1));
    //   j = 1/2 ... ny-1/2
    for (I j = 0; j < ny; ++j) {
      auto y = y0 + dy * (2*j + 1) / 2;
      g += 8 * func(x, y);
    }
    //   j = 1 ... ny-1
    for (I j = 1; j < ny; ++j) {
      auto y = y0 + dy * j;
      g += 4 * func(x, y);
    }
  }
  // i == nx
  {
    auto x = x1;
    //   j == 0 || j == ny
    g += func(x, y0) + func(x, y1);
    //   j = 1/2 ... ny-1/2
    for (I j = 0; j < ny; ++j) {
      auto y = y0 + dy * (2*j + 1) / 2;
      g += 4 * func(x, y);
    }
    //   j = 1 ... ny-1
    for (I j = 1; j < ny; ++j) {
      auto y = y0 + dy * j;
      g += 2 * func(x, y);
    }
  }
  g *= dx * dy / 36;
  return g;
}

} // end namespace standards
