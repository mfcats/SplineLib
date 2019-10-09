/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "gmock/gmock.h"

#include "src/util/vector_utils.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

using namespace splinelib::src;

class Vectors : public Test {
 public:
  Vectors() : int_vectorA({2, 3, 4, 5}),
              int_vectorB({3, 2, 1, 0}),
              double_vectorA({4.5, 6.7, 8.9}),
              double_vectorB({3.2, 4.3, 5.4}) {}

 protected:
  std::vector<int> int_vectorA;
  std::vector<int> int_vectorB;
  std::vector<double> double_vectorA;
  std::vector<double> double_vectorB;
};

TEST_F(Vectors, CanComputeTwoNorm) {  // NOLINT
  ASSERT_THAT(util::VectorUtils<int>::ComputeTwoNorm(int_vectorB), DoubleEq(sqrt(14)));
  ASSERT_THAT(util::VectorUtils<double>::ComputeTwoNorm(double_vectorA), DoubleEq(sqrt(144.35)));
}

TEST_F(Vectors, CanComputeDifference) {  // NOLINT
  ASSERT_THAT(util::VectorUtils<int>::ComputeDifference(int_vectorA, int_vectorB)[0], -1);
  ASSERT_THAT(util::VectorUtils<int>::ComputeDifference(int_vectorA, int_vectorB)[1], 1);
  ASSERT_THAT(util::VectorUtils<int>::ComputeDifference(int_vectorA, int_vectorB)[2], 3);
  ASSERT_THAT(util::VectorUtils<int>::ComputeDifference(int_vectorA, int_vectorB)[3], 5);
  ASSERT_THAT(util::VectorUtils<double>::ComputeDifference(double_vectorA, double_vectorB)[0], DoubleEq(1.3));
  ASSERT_THAT(util::VectorUtils<double>::ComputeDifference(double_vectorA, double_vectorB)[1], DoubleEq(2.4));
  ASSERT_THAT(util::VectorUtils<double>::ComputeDifference(double_vectorA, double_vectorB)[2], DoubleEq(3.5));
}

TEST_F(Vectors, CanComputeDistance) {  // NOLINT
  ASSERT_THAT(util::VectorUtils<int>::ComputeDistance(int_vectorA, int_vectorB), DoubleEq(6));
  ASSERT_THAT(util::VectorUtils<double>::ComputeDistance(double_vectorA, double_vectorB), DoubleEq(sqrt(19.7)));
}

TEST_F(Vectors, CanComputeScalarProduct) {  // NOLINT
  ASSERT_THAT(util::VectorUtils<int>::ComputeScalarProduct(int_vectorA, int_vectorB), 16);
  ASSERT_THAT(util::VectorUtils<double>::ComputeScalarProduct(double_vectorA, double_vectorB), DoubleEq(91.27));
}

TEST_F(Vectors, CanScaleVector) {  // NOLINT
  ASSERT_THAT(util::VectorUtils<int>::ScaleVector(int_vectorA, 2)[0], 4);
  ASSERT_THAT(util::VectorUtils<int>::ScaleVector(int_vectorA, 2)[1], 6);
  ASSERT_THAT(util::VectorUtils<int>::ScaleVector(int_vectorA, 2)[2], 8);
  ASSERT_THAT(util::VectorUtils<int>::ScaleVector(int_vectorA, 2)[3], 10);
  ASSERT_THAT(util::VectorUtils<double>::ScaleVector(double_vectorA, 0.5)[0], DoubleEq(2.25));
  ASSERT_THAT(util::VectorUtils<double>::ScaleVector(double_vectorA, 0.5)[1], DoubleEq(3.35));
  ASSERT_THAT(util::VectorUtils<double>::ScaleVector(double_vectorA, 0.5)[2], DoubleEq(4.45));
}

TEST_F(Vectors, CanComputeCrossProduct) {  // NOLINT
  ASSERT_THAT(util::VectorUtils<double>::CrossProduct(double_vectorA, double_vectorB)[0], DoubleNear(-2.09, 1e-10));
  ASSERT_THAT(util::VectorUtils<double>::CrossProduct(double_vectorA, double_vectorB)[1], DoubleEq(4.18));
  ASSERT_THAT(util::VectorUtils<double>::CrossProduct(double_vectorA, double_vectorB)[2], DoubleNear(-2.09, 1e-10));
}

TEST_F(Vectors, CanBeFiltered) {  // NOLINT
  std::vector<int> filtered = util::VectorUtils<int>::FilterVector(int_vectorA, {1, 2});
  ASSERT_THAT(filtered.size(), 2);
  ASSERT_THAT(filtered[0], 3);
  ASSERT_THAT(filtered[1], 4);
  ASSERT_THROW(util::VectorUtils<double>::FilterVector(double_vectorA, {2, 3}), std::runtime_error);
}
