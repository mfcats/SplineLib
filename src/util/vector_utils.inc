/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

// NOLINT(build/header_guard)

#include <algorithm>
#include <functional>
#include <vector>

template<typename TYPE>
std::vector<TYPE> ScaleVector(std::vector<TYPE> const &vectorA, TYPE factor) {
  std::vector<TYPE> result(vectorA);
  std::transform(result.begin(), result.end(), result.begin(),
                 std::bind(std::multiplies<TYPE>(), factor, std::placeholders::_1));
  return result;
}

template<typename TYPE>
TYPE ComputeScalarProduct(std::vector<TYPE> const &vectorA, std::vector<TYPE> const &vectorB) {
  return std::inner_product(vectorA.begin(), vectorA.end(), vectorB.begin(), 0.0);
}

template<typename TYPE>
std::vector<TYPE> ComputeSum(std::vector<TYPE> const &vectorA, std::vector<TYPE> const &vectorB) {
  std::vector<TYPE> sum(vectorB);
  std::transform(vectorA.begin(), vectorA.end(), sum.begin(), sum.begin(), std::plus<TYPE>());
  return sum;
}

template<typename TYPE>
std::vector<TYPE> ComputeDifference(std::vector<TYPE> const &vectorA, std::vector<TYPE> const &vectorB) {
  std::vector<TYPE> difference(vectorB);
  std::transform(vectorA.begin(), vectorA.end(), difference.begin(), difference.begin(), std::minus<TYPE>());
  return difference;
}

template<typename TYPE>
double ComputeTwoNorm(std::vector<TYPE> const &vectorA) {
  return sqrt(ComputeScalarProduct(vectorA, vectorA));
}

template<typename TYPE>
double ComputeDistance(std::vector<TYPE> const &vectorA, std::vector<TYPE> const &vectorB) {
  return ComputeTwoNorm<TYPE>(ComputeDifference<TYPE>(vectorA, vectorB));
}

template<typename TYPE>
std::vector<TYPE> ComputeCrossProduct(std::vector<TYPE> const &a, std::vector<TYPE> const &b) {
  std::vector<TYPE> result(a.size());
  result[0] = a[1] * b[2] - a[2] * b[1];
  result[1] = a[2] * b[0] - a[0] * b[2];
  result[2] = a[0] * b[1] - a[1] * b[0];
  return result;
}

template<typename TYPE>
std::vector<TYPE> GetEntriesAtIndices(std::vector<TYPE> const &input, std::vector<int> const &indices) {
  std::vector<TYPE> output;
  for (const auto &index : indices) {
    if (index < static_cast<int>(input.size())) {
      output.emplace_back(input[index]);
    } else {
      throw std::runtime_error("The vector index is too high to be filtered from the input vector.");
    }
  }
  return output;
}
