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
#include <cmath>
#include <functional>
#include <numeric>
#include <vector>

namespace splinelib::src::util::vector_utils {
template<typename TYPE>
double ComputeTwoNorm(std::vector<TYPE> vectorA) {
  std::transform(vectorA.begin(), vectorA.end(), vectorA.begin(), vectorA.begin(), std::multiplies<TYPE>());
  double sum = 0;
  for (TYPE i : vectorA) {
    sum += i;
  }
  return sqrt(sum);
}

template<typename TYPE>
std::vector<TYPE> ComputeDifference(std::vector<TYPE> vectorA, std::vector<TYPE> vectorB) {
  std::transform(vectorA.begin(), vectorA.end(), vectorB.begin(), vectorB.begin(), std::minus<TYPE>());
  return vectorB;
}

template<typename TYPE>
double ComputeDistance(std::vector<TYPE> vectorA, std::vector<TYPE> vectorB) {
  return ComputeTwoNorm<TYPE>(ComputeDifference<TYPE>(vectorA, vectorB));
}

template<typename TYPE>
TYPE ComputeScalarProduct(std::vector<TYPE> vectorA, std::vector<TYPE> vectorB) {
  std::transform(vectorA.begin(), vectorA.end(), vectorB.begin(), vectorB.begin(), std::multiplies<TYPE>());
  TYPE sum = 0;
  for (TYPE i : vectorB) {
    sum += i;
  }
  return sum;
}

template<typename TYPE>
std::vector<TYPE> ScaleVector(std::vector<TYPE> vectorA, TYPE factor) {
  std::transform(vectorA.begin(), vectorA.end(), vectorA.begin(),
                 std::bind(std::multiplies<TYPE>(), factor, std::placeholders::_1));
  return vectorA;
}

template<typename TYPE>
std::vector<TYPE> CrossProduct(std::vector<TYPE> const &a, std::vector<TYPE> const &b) {
  std::vector<TYPE> r(a.size());
  r[0] = a[1] * b[2] - a[2] * b[1];
  r[1] = a[2] * b[0] - a[0] * b[2];
  r[2] = a[0] * b[1] - a[1] * b[0];
  return r;
}

template<typename TYPE>
std::vector<TYPE> FilterVector(const std::vector<TYPE> &input, const std::vector<int> &positions) {
  std::vector<TYPE> output;
  for (const auto &pos : positions) {
    if (pos < static_cast<int>(input.size())) {
      output.emplace_back(input[pos]);
    } else {
      throw std::runtime_error("The vector index is too high to be filtered from the input vector.");
    }
  }
  return output;
}

#include "src/util/vector_utils.inc"
}  // namespace splinelib::src::util::vector_utils

#endif  // SRC_UTIL_VECTOR_UTILS_H_
