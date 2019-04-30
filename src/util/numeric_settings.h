/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_UTIL_NUMERIC_SETTINGS_H_
#define SRC_UTIL_NUMERIC_SETTINGS_H_

#include <cmath>
#include <limits>

#include "named_type.h"

namespace util {
template<typename T>
class NumericSettings {
 public:
  constexpr static T kEpsilon() {
    return kEpsilonFactor_ * std::numeric_limits<T>::epsilon();
  }

  constexpr static bool AreEqual(const T &a, const T &b, T tolerance = kEpsilon()) {
    return std::abs(a - b) < tolerance;
  }

  NumericSettings(const NumericSettings &numericSettings) = delete;

  void operator=(const NumericSettings &numericSettings) = delete;

 private:
  NumericSettings() = default;

  constexpr static double kEpsilonFactor_ = 10.0;
};
}  // namespace util

#endif  // SRC_UTIL_NUMERIC_SETTINGS_H_
