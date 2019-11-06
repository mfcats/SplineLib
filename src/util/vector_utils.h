/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_UTIL_VECTOR_UTILS_H_
#define SRC_UTIL_VECTOR_UTILS_H_

#include <cmath>
#include <numeric>
#include <vector>

namespace splinelib::src::util::vector_utils {
template<typename TYPE>
std::vector<TYPE> ScaleVector(std::vector<TYPE> const &vectorA, TYPE factor);

template<typename TYPE>
std::vector<TYPE> ComputeSum(std::vector<TYPE> const &vectorA, std::vector<TYPE> const &vectorB);
template<typename TYPE>
std::vector<TYPE> ComputeDifference(std::vector<TYPE> const &vectorA, std::vector<TYPE> const &vectorB);

template<typename TYPE>
TYPE ComputeScalarProduct(std::vector<TYPE> const &vectorA, std::vector<TYPE> const &vectorB);
template<typename TYPE>
double ComputeTwoNorm(std::vector<TYPE> const &vectorA);
template<typename TYPE>
double ComputeDistance(std::vector<TYPE> const &vectorA, std::vector<TYPE> const &vectorB);
template<typename TYPE>
std::vector<TYPE> ComputeCrossProduct(std::vector<TYPE> const &a, std::vector<TYPE> const &b);

template<typename TYPE>
std::vector<TYPE> GetEntriesAtIndices(std::vector<TYPE> const &input, std::vector<int> const &indices);

#include "src/util/vector_utils.inc"
}  // namespace splinelib::src::util::vector_utils

#endif  // SRC_UTIL_VECTOR_UTILS_H_
