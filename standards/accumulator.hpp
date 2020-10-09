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

// Ref: http://modelingwithdata.org/pdfs/moments.pdf
// Ref: standards/doc/stat/fourth-moment.pdf
// Ref: wrn:2017-01-20, 2017-02-07, 2017-02-09

#pragma once

#include <cmath>
#include <string>
#include "power.hpp"
#include "moment.hpp"

namespace standards {

class accumulator : public moment<accumulator> {
private:
  typedef moment<accumulator> super_type;
public:
  accumulator(std::string const& name = "") : super_type(*this), name_(name) { reset(); }
  accumulator(accumulator const& x) : super_type(*this), name_(x.name_),
    count_(x.count_), sum1_(x.sum1_), sum2_(x.sum2_), sum3_(x.sum3_), sum4_(x.sum4_) {}
  void reset() {
    count_ = 0;
    sum1_ = sum2_ = sum3_ = sum4_ = 0;
  }
  double operator<<(double v) {
    ++count_;
    sum1_ += v;
    sum2_ += p2(v);
    sum3_ += p3(v);
    sum4_ += p4(v);
    return v;
  }

  void set_name(std::string const& name) { name_ = name; }
  std::string name() const { return name_; }
  unsigned long count() const { return count_; }
  double moment1() const { return count() ? (sum1_ / count()) : 0; }
  double moment2() const { return count() ? (sum2_ / count()) : 0; }
  double moment3() const { return count() ? (sum3_ / count()) : 0; }
  double moment4() const { return count() ? (sum4_ / count()) : 0; }

  double central_moment1() const { return 0; }
  double central_moment2() const {
    double n = count();
    return (n > 1) ? (n/(n-1)) * super_type::central_moment2() : 0;
  }
  double central_moment3() const {
    double n = count();
    return (n > 2) ? (n*n/((n-1)*(n-2))) * super_type::central_moment3() : 0;
  }
  double central_moment4() const {
    double n = count();
    return (n > 3) ?
      n/((n-1)*(n-2)*(n-3)) *
      ((p2(n)-2*n+3)*super_type::central_moment4() -
       (6*n-9)*(p2(moment2()) - 2 * p2(moment1()) * moment2() + p4(moment1()))) : 0;
  }
  double cumulant4() const {
    double n = count();
    return (n > 3) ?
      p2(n)/((n-1)*(n-2)*(n-3)) *
      ((n+1)*super_type::central_moment4() -
       3*(n-1)*(p2(moment2()) - 2 * p2(moment1()) * moment2() + p4(moment1()))) : 0;
  }
  double average() const { return super_type::mean(); }
  double error() const { return count() ? std::sqrt(super_type::variance() / count()) : 0; }
private:
  std::string name_;
  unsigned long count_;
  double sum1_, sum2_, sum3_, sum4_;
};

inline std::string format(accumulator const& accum) {
  return accum.name() + " = "
    + std::to_string(accum.average()) + " +- "
    + std::to_string(accum.error());
}

inline std::ostream& operator<<(std::ostream& os, accumulator const& accum) {
  os << accum.name() << " = " << accum.average() << " +- " << accum.error();
  return os;
}

} // end namespace standards
