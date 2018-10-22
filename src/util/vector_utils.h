/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_UTIL_VECTOR_UTILS_H_
#define SRC_UTIL_VECTOR_UTILS_H_

#include <algorithm>
#include <functional>
#include <vector>

namespace util {
template<typename T>
class VectorUtils {
 public:
  static T ComputeTwoNorm(std::vector<T> vectorA) {
    std::transform(vectorA.begin(), vectorA.end(), vectorA.begin(), vectorA.begin(), std::multiplies<T>());
    T sum = 0;
    for (T i : vectorA) {
      sum += i;
    }
    return sqrt(sum);
  }

  static std::vector<T> ComputeDifference(std::vector<T> vectorA, std::vector<T> vectorB) {
    std::transform(vectorA.begin(), vectorA.end(), vectorB.begin(), vectorB.begin(), std::minus<T>());
    return vectorB;
  }

  static T ComputeScalarProduct(std::vector<T> vectorA, std::vector<T> vectorB) {
    std::transform(vectorA.begin(), vectorA.end(), vectorB.begin(), vectorB.begin(), std::multiplies<T>());
    return std::accumulate(vectorB.begin(), vectorB.end(), 0);
  }

  static std::vector<T> ScaleVector(std::vector<T> vectorA, T factor) {
    std::transform(vectorA.begin(), vectorA.end(), vectorA.begin(), std::bind1st(std::multiplies<T>(), factor));
    return vectorA;
  }

  static std::vector<T> CrossProduct(std::vector<T> const &a, std::vector<T> const &b) {
    std::vector<T> r(a.size());
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
    return r;
  }
};
}  // namespace util

#endif  // SRC_UTIL_VECTOR_UTILS_H_
