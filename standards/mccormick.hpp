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

// McCormick function, which takes minimum -1.91322 at (-0.547198, -1.5472)
  
class mccormick {
public:
  double operator()(std::vector<double> const& x) const {
    return std::sin(x[0]+x[1]) + (x[0]-x[1]) * (x[0]-x[1])
      - 1.5 * x[0] + 2.5 * x[1] + 1;
  }
};

} // end namespace standards
