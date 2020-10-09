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

namespace standards {

class golden_section {
public:
  // one (or more) local minimum should be included in (x0,x2)
  template<class FUNC>
  int find_minimum(FUNC& f, double x0, double x2,
    double prec = std::sqrt(2 * std::numeric_limits<double>::epsilon())) {
    if (x2 < x0) std::swap(x0, x2);
    double y0 = f(x0);
    double y2 = f(x2);
    int counter = 2;
    // initial enclosure
    double dx = x2 - x0;
    double x1 = x0 + 0.5 * dx;
    double y1;
    while (true) {
      ++counter;
      y1 = f(x1);
      if (y1 < y0 && y1 < y2) break;
      x1 += dx;
      if (x1 > x2) {
        dx *= 0.5;
        x1 = x0 + 0.5 * dx;
      }
      if (counter > 1024) return -1;
    }
    // start golden section
    double r = (3-std::sqrt(5))/2; // golden ratio
    double prec2 = prec * prec;
    while (((x2-x0) / std::max(std::abs(x0), std::abs(x2))) > prec ||
           ((y0-y1) / std::max(std::abs(y1), std::abs(y0))) > prec2 ||
           ((y2-y1) / std::max(std::abs(y1), std::abs(y2))) > prec2) {
      ++counter;
      if ((x2-x1) > (x1-x0)) {
        double xn = (x2-x1) * r + x1;
        double yn = f(xn);
        if (y1 > yn) {
          x0 = x1; y0 = y1;
          x1 = xn; y1 = yn;
        } else {
          x2 = xn; y2 = yn;
        }
      } else{
        double xn = (x0-x1) * r + x1;
        double yn = f(xn);
        if (y1 > yn) {
          x2 = x1; y2 = y1;
          x1 = xn; y1 = yn;
        } else {
          x0 = xn; y0 = yn;
        }
      }
    }
    minarg_ = x1;
    minval_ = y1;
    return counter;
  }
  double minarg() const { return minarg_; }
  double minval() const { return minval_; }
private:
  double minarg_, minval_;
};

} // end namespace standards
