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
#include <string>
#include "moment.hpp"

namespace standards {

class bernoulli_distribution : public moment<bernoulli_distribution> {
private:
  typedef moment<bernoulli_distribution> super_type;
public:
  bernoulli_distribution(double p) : super_type(*this), p_(p) {
    if (p < 0 || p >1)
      throw std::invalid_argument("standards::bernoulli_distribution");
  }
  std::string name() const {
    return "Bernoulli Distribution: Be(" + std::to_string(p_) + ")";
  }
  double moment1() const { return p_; }
  double moment2() const { return p_; }
  double moment3() const { return p_; }
  double moment4() const { return p_; }
private:
  double p_;
};

} // end namespace standards
