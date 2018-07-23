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
    return sqrt(std::accumulate(vectorA.begin(), vectorA.end(), 0));
  }

  static std::vector<double> ComputeDifference(std::vector<double> vectorA, std::vector<double> vectorB) {
    std::transform(vectorA.begin(), vectorA.end(), vectorB.begin(), vectorB.begin(), std::minus<double>());
    return vectorB;
  }

  static double ComputeScalarProduct(std::vector<double> vectorA, std::vector<double> vectorB) {
    std::transform(vectorA.begin(), vectorA.end(), vectorB.begin(), vectorB.begin(), std::multiplies<double>());
    return std::accumulate(vectorB.begin(), vectorB.end(), 0);
  }
};
}  // namespace util

#endif  // SRC_UTIL_VECTOR_UTILS_H_
