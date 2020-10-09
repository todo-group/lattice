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
#include "power.hpp"

namespace standards {

template<typename DIST>
class moment {
private:
  typedef DIST dist_type;
public:
  moment(dist_type const& dist) : dist_(dist) {}

  double central_moment1() const { return 0; }
  double central_moment2() const { return dist_.moment2() - p2(dist_.moment1()); }
  double central_moment3() const {
    return dist_.moment3() - 3 * dist_.moment2() * dist_.moment1() + 2 * p3(dist_.moment1());
  }
  double central_moment4() const {
    return dist_.moment4() - 4 * dist_.moment3() * dist_.moment1() 
      + 6 * dist_.moment2() * p2(dist_.moment1()) - 3 * p4(dist_.moment1());
  }

  double cumulant1() const { return dist_.moment1(); }
  double cumulant2() const { return dist_.central_moment2(); }
  double cumulant3() const { return dist_.central_moment3(); }
  double cumulant4() const {
    return dist_.central_moment4() - 3 * p2(dist_.central_moment2());
  }

  double mean() const { return dist_.cumulant1(); }
  double variance() const { return dist_.cumulant2(); }
  double standard_deviation() const { return std::sqrt(variance()); }
  double skewness() const {
    return (standard_deviation() > 0) ? dist_.central_moment3() / p3(standard_deviation()) : 0;
  }
  double kurtosis() const {
    return (standard_deviation() > 0) ? dist_.central_moment4() / p4(standard_deviation()) : 0;
  }
  double kurtosis_excess() const { return kurtosis() - 3; }
private:
  dist_type const& dist_;
};

} // end namespace standards
