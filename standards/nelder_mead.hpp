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

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace standards {

class nelder_mead {
public:
  template<class FUNC>
  int find_minimum(FUNC& f, std::vector<double> x0, double shift = 0.1,
    double prec = 2 * std::numeric_limits<double>::epsilon(),
    int max_iteration = 1024) {
    int n = x0.size();
    x.resize(n+1);
    for (int i = 0; i < n+1; ++i) x[i].resize(n);
    y.resize(n+1);
    xn.resize(n); // newly proposed point
    xg.resize(n); // cener of mass of x[0]...x[n-1]
    // initial condition
    x[0] = x0;
    y[0] = f(x[0]);
    for (int i = 0; i < n; ++i) {
      x[i+1] = x0;
      x[i+1][i] += shift;
      y[i+1] = f(x[i+1]);
    }
    int counter = n+1;
    while (true) {
      // set y[n] as the largest, y[n-1] as the second largest, and y[0] the smallest
      for (int i = 0; i < n; ++i) {
        if (y[i] > y[n]) {
          std::swap(x[i], x[n]);
          std::swap(y[i], y[n]);
        }
      }
      for (int i = 1; i < n; ++i) {
        if (y[i] < y[0]) {
          std::swap(x[i], x[0]);
          std::swap(y[i], y[0]);
        }
      }
      for (int i = 1; i < n-1; ++i) {
        if (y[i] > y[n-1]) {
          std::swap(x[i], x[n-1]);
          std::swap(y[i], y[n-1]);
        }
      }
      // center of mass except for x[n]
      for (int j = 0; j < n; ++j) xg[j] = 0.0;
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          xg[j] += x[i][j];
        }
      }
      for (int j = 0; j < n; ++j) xg[j] /= n;
      // refrection
      for (int j = 0; j < n; ++j) xn[j] = xg[j] - (x[n][j] - xg[j]);
      ++counter;
      double yn = f(xn);
      if (yn < y[n]){
        for (int j = 0; j < n; ++j) x[n][j] = xn[j];
        y[n] = yn;
      }
      if (yn < y[0]) {
        // reduction and expansion
        for (int j = 0; j < n; ++j) xn[j] = xg[j] + 2 * (x[n][j] - xg[j]);
        ++counter;
        yn = f(xn);
        if (yn < y[0]) {
          for (int j = 0; j < n; ++j) x[n][j] = xn[j];
          y[n] = yn;
        }
      } else if (y[n] > y[n-1]) {
        // contraction
        for (int j = 0; j < n; ++j) xn[j] = xg[j] + 0.5 * (x[n][j] - xg[j]);
        ++counter;
        yn = f(xn);
        if (yn < y[n]) {
          for (int j = 0; j < n; ++j) x[n][j] = xn[j];
          y[n] = yn;
        } else {
          // contraction in all directions
          for (int i = 1; i <= n; ++i) {
            ++counter;
            for (int j = 0; j < n; ++j) x[i][j] = x[0][j] + 0.5 * (x[i][j] - x[0][j]);
            y[i] = f(x[i]);
          }
        }
      }
      // convergence check
      double diff2 = 0.0;
      for (int j = 0; j < n; ++j) diff2 += (x[n][j] - x[0][j]) * (x[n][j] - x[0][j]);
      if (diff2 < prec || std::abs(y[n] - y[0]) < prec * std::max(std::abs(y[n]), std::abs(y[0]))) {
        return counter;
      }
      if (counter > max_iteration) return -2;
    }
  }
  std::vector<double> const& minarg() const { return x[0]; }
  double minval() const { return y[0]; }
private:
  std::vector<std::vector<double> > x;
  std::vector<double> y;
  std::vector<double> xn, xg;
};

} // end namespace standards
