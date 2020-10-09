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

namespace standards {

class rastrigin {
public:
  rastrigin(unsigned int n) : n_(n) {}
  double operator()(std::vector<double> const& x) const {
    const double a = 10;
    double f = a * n_;
    for (unsigned int i = 0; i < n_; ++i) {
      f += x[i] * x[i] - a * std::cos(2 * M_PI * x[i]);
    }
    return f;
  }
private:
  double n_;
};

} // end namespace standards
