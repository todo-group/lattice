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

#include <chrono>

namespace standards {

class timer {
public:
  timer() : start_(std::chrono::system_clock::now()) {}
  void reset() {
    start_ = std::chrono::system_clock::now();
  }
  double elapsed() const {
    auto end = std::chrono::system_clock::now();
    return 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(end - start_).count();
  }
private:
  std::chrono::system_clock::time_point start_;
};

} // end namespace standards
