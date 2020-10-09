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

#include <stdexcept>
#include <string>
#include "power.hpp"
#include "moment.hpp"

namespace standards {

class normal_distribution : public moment<normal_distribution> {
private:
  typedef moment<normal_distribution> super_type;
public:
  normal_distribution(double mu = 0, double sigma = 1) : super_type(*this), mu_(mu), sigma_(sigma) {
    if (sigma_ <= 0)
      throw std::invalid_argument("standards::normal_distribution");
  }
  std::string name() const {
    return "Normal Distribution: N(" + std::to_string(mu_) + ","
      + std::to_string(sigma_) + ")";
  }
  double moment1() const { return mu_; }
  double moment2() const { return p2(mu_) + p2(sigma_); }
  double moment3() const { return p3(mu_) + 3 * mu_ * p2(sigma_); }
  double moment4() const { return p4(mu_) + 6 * p2(mu_) * p2(sigma_) + 3 * p4(sigma_); }
private:
  double mu_, sigma_;
};

} // end namespace standards
