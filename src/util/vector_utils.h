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
  static double ComputeTwoNorm(std::vector<double> vectorA) {
    std::transform(vectorA.begin(), vectorA.end(), vectorA.begin(), vectorA.begin(), std::multiplies<double>());
    auto r = vectorA[0];
    auto s = vectorA[1];
    auto t = vectorA[2];
    auto v = r + s;
    auto w = v + t;
    double u = std::accumulate(vectorA.begin(), vectorA.end(), 0);
    return sqrt(w);
  }

  static std::vector<double> ComputeDifference(std::vector<double> vectorA, std::vector<double> vectorB) {
    std::transform(vectorA.begin(), vectorA.end(), vectorB.begin(), vectorB.begin(), std::minus<double>());
    return vectorB;
  }

  static double ComputeScalarProduct(std::vector<double> vectorA, std::vector<double> vectorB) {
    std::transform(vectorA.begin(), vectorA.end(), vectorB.begin(), vectorB.begin(), std::multiplies<double>());
    return std::accumulate(vectorB.begin(), vectorB.end(), 0);
  }

  static std::vector<double> ScaleVector(std::vector<double> vectorA, double factor) {
    std::transform(vectorA.begin(), vectorA.end(), vectorA.begin(), std::bind1st(std::multiplies<T>(), factor));
    return vectorA;
  }

  static std::vector<double> CrossProduct(std::vector<double> const &a, std::vector<double> const &b) {
    std::vector<double> r(a.size());
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
    return r;
  }
};
}  // namespace util

#endif  // SRC_UTIL_VECTOR_UTILS_H_
