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

class bisection {
public:
  // zero(s) should be included in (x0,x1)
  template<class FUNC>
  int find_zero(FUNC& f, double x0, double x1,
    double prec = 2 * std::numeric_limits<double>::epsilon()) {
    if (x0 > x1) std::swap(x0, x1);
    double y0 = f(x0);
    int counter = 1;
    if (std::abs(y0) < prec) {
      zero_ = x0;
      return counter;
    }
    ++counter;
    double y1 = f(x1);
    if (std::abs(y1) < prec) {
      zero_ = x1;
      return counter;
    }
    if (y0 * y1 > 0) return -1;
    // start bisection
    while (((x1-x0) / std::max(std::abs(x0), std::abs(x1))) > prec &&
           std::max(std::abs(y1), std::abs(y0)) > prec) {
      ++counter;
      double xn = (x0 + x1) / 2;
      double yn = f(xn);
      if (y0 * yn > 0) {
        x0 = xn; y0 = yn;
      } else{
        x1 = xn; y1 = yn;
      }
    }
    zero_ = x0;
    return counter;
  }
  double zero() const { return zero_; }
private:
  double zero_;
};

} // end namespace standards
