/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/util/vector_utils.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

using namespace splinelib::src;

class TwoIntegerVectorsForVectorUtils : public Test {
 public:
  TwoIntegerVectorsForVectorUtils() : int_vector_a_({2, 3, 4, 5}), int_vector_b_({3, 2, 1, 0}) {}

 protected:
  std::vector<int> int_vector_a_;
  std::vector<int> int_vector_b_;
};

TEST_F(TwoIntegerVectorsForVectorUtils, CanBeScaledWithFactor2) {  // NOLINT
  ASSERT_THAT(util::vector_utils::ScaleVector(int_vector_a_, 2), std::vector({4, 6, 8, 10}));
  ASSERT_THAT(util::vector_utils::ScaleVector(int_vector_b_, 2), std::vector({6, 4, 2, 0}));
}

TEST_F(TwoIntegerVectorsForVectorUtils, CanBeAdded) {  // NOLINT
  ASSERT_THAT(util::vector_utils::ComputeSum(int_vector_a_, int_vector_b_), std::vector({5, 5, 5, 5}));
  ASSERT_THAT(util::vector_utils::ComputeSum(int_vector_b_, int_vector_a_), std::vector({5, 5, 5, 5}));
}

TEST_F(TwoIntegerVectorsForVectorUtils, CanBeSubtracted) {  // NOLINT
  ASSERT_THAT(util::vector_utils::ComputeDifference(int_vector_a_, int_vector_b_), std::vector({-1, 1, 3, 5}));
  ASSERT_THAT(util::vector_utils::ComputeDifference(int_vector_b_, int_vector_a_), std::vector({1, -1, -3, -5}));
}

TEST_F(TwoIntegerVectorsForVectorUtils, CanComputeScalarProductOf16) {  // NOLINT
  ASSERT_THAT(util::vector_utils::ComputeScalarProduct(int_vector_a_, int_vector_b_), 16);
  ASSERT_THAT(util::vector_utils::ComputeScalarProduct(int_vector_b_, int_vector_a_), 16);
}

TEST_F(TwoIntegerVectorsForVectorUtils, CanComputeTwoNormOfsqrt54Andsqrt14) {  // NOLINT
  ASSERT_THAT(util::vector_utils::ComputeTwoNorm(int_vector_a_), DoubleEq(sqrt(54)));
  ASSERT_THAT(util::vector_utils::ComputeTwoNorm(int_vector_b_), DoubleEq(sqrt(14)));
}

TEST_F(TwoIntegerVectorsForVectorUtils, CanComputeDistanceOf6) {  // NOLINT
  ASSERT_THAT(util::vector_utils::ComputeDistance(int_vector_a_, int_vector_b_), DoubleEq(6));
  ASSERT_THAT(util::vector_utils::ComputeDistance(int_vector_b_, int_vector_a_), DoubleEq(6));
}

TEST_F(TwoIntegerVectorsForVectorUtils, CanComputeCrossProduct) {  // NOLINT
  std::vector<int> int_vector_a_with_first_three_components{int_vector_a_.begin(), int_vector_a_.begin() + 3};
  std::vector<int> int_vector_b_with_first_three_components{int_vector_b_.begin(), int_vector_b_.begin() + 3};
  ASSERT_THAT(util::vector_utils::ComputeCrossProduct(
      int_vector_a_with_first_three_components, int_vector_b_with_first_three_components), std::vector({-5, 10, -5}));
  ASSERT_THAT(util::vector_utils::ComputeCrossProduct(
      int_vector_b_with_first_three_components, int_vector_a_with_first_three_components), std::vector({5, -10, 5}));
}

TEST_F(TwoIntegerVectorsForVectorUtils, CanBeFiltered) {  // NOLINT
  std::vector<int> filtered = util::vector_utils::GetEntriesAtIndices(int_vector_a_, {1, 2});
  ASSERT_THAT(filtered.size(), 2);
  ASSERT_THAT(filtered[0], 3);
  ASSERT_THAT(filtered[1], 4);
  ASSERT_THROW(util::vector_utils::GetEntriesAtIndices(int_vector_b_, {2, 4}), std::invalid_argument);
}

class TwoDoubleVectorsForVectorUtils : public Test {
 public:
  TwoDoubleVectorsForVectorUtils() : double_vector_a_({4.5, 6.7, 8.9}), double_vector_b_({3.2, 4.3, 5.4}) {}

 protected:
  std::vector<double> double_vector_a_;
  std::vector<double> double_vector_b_;
};

TEST_F(TwoDoubleVectorsForVectorUtils, CanBeScaledWithFactorMinus0_5) {  // NOLINT
  ASSERT_THAT(util::vector_utils::ScaleVector(double_vector_a_, -0.5), std::vector({-2.25, -3.35, -4.45}));
  ASSERT_THAT(util::vector_utils::ScaleVector(double_vector_b_, -0.5), std::vector({-1.6, -2.15, -2.7}));
}

TEST_F(TwoDoubleVectorsForVectorUtils, CanBeAdded) {  // NOLINT
  ASSERT_THAT(util::vector_utils::ComputeSum(double_vector_a_, double_vector_b_), std::vector({7.7, 11.0, 14.3}));
  ASSERT_THAT(util::vector_utils::ComputeSum(double_vector_b_, double_vector_a_), std::vector({7.7, 11.0, 14.3}));
}

TEST_F(TwoDoubleVectorsForVectorUtils, CanBeSubtracted) {  // NOLINT
  ASSERT_THAT(util::vector_utils::ComputeDifference(double_vector_a_, double_vector_b_)[0], DoubleNear(1.3, 1e-10));
  ASSERT_THAT(util::vector_utils::ComputeDifference(double_vector_a_, double_vector_b_)[1], DoubleNear(2.4, 1e-10));
  ASSERT_THAT(util::vector_utils::ComputeDifference(double_vector_a_, double_vector_b_)[2], DoubleEq(3.5));
  ASSERT_THAT(util::vector_utils::ComputeDifference(double_vector_b_, double_vector_a_)[0], DoubleNear(-1.3, 1e-10));
  ASSERT_THAT(util::vector_utils::ComputeDifference(double_vector_b_, double_vector_a_)[1], DoubleNear(-2.4, 1e-10));
  ASSERT_THAT(util::vector_utils::ComputeDifference(double_vector_b_, double_vector_a_)[2], DoubleEq(-3.5));
}

TEST_F(TwoDoubleVectorsForVectorUtils, CanComputeScalarProductOf91_27) {  // NOLINT
  ASSERT_THAT(util::vector_utils::ComputeScalarProduct(double_vector_a_, double_vector_b_), DoubleEq(91.27));
  ASSERT_THAT(util::vector_utils::ComputeScalarProduct(double_vector_b_, double_vector_a_), DoubleEq(91.27));
}

TEST_F(TwoDoubleVectorsForVectorUtils, CanComputeTwoNormOfsqrt144_35Andsqrt57_89) {  // NOLINT
  ASSERT_THAT(util::vector_utils::ComputeTwoNorm(double_vector_a_), DoubleEq(sqrt(144.35)));
  ASSERT_THAT(util::vector_utils::ComputeTwoNorm(double_vector_b_), DoubleEq(sqrt(57.89)));
}

TEST_F(TwoDoubleVectorsForVectorUtils, CanComputeDistanceOf19_7) {  // NOLINT
  ASSERT_THAT(util::vector_utils::ComputeDistance(double_vector_a_, double_vector_b_), DoubleEq(sqrt(19.7)));
  ASSERT_THAT(util::vector_utils::ComputeDistance(double_vector_b_, double_vector_a_), DoubleEq(sqrt(19.7)));
}

TEST_F(TwoDoubleVectorsForVectorUtils, CanComputeCrossProduct) {  // NOLINT
  ASSERT_THAT(util::vector_utils::ComputeCrossProduct(double_vector_a_, double_vector_b_)[0], DoubleNear(-2.09, 1e-10));
  ASSERT_THAT(util::vector_utils::ComputeCrossProduct(double_vector_a_, double_vector_b_)[1], DoubleEq(4.18));
  ASSERT_THAT(util::vector_utils::ComputeCrossProduct(double_vector_a_, double_vector_b_)[2], DoubleNear(-2.09, 1e-10));
  ASSERT_THAT(util::vector_utils::ComputeCrossProduct(double_vector_b_, double_vector_a_)[0], DoubleNear(2.09, 1e-10));
  ASSERT_THAT(util::vector_utils::ComputeCrossProduct(double_vector_b_, double_vector_a_)[1], DoubleEq(-4.18));
  ASSERT_THAT(util::vector_utils::ComputeCrossProduct(double_vector_b_, double_vector_a_)[2], DoubleNear(2.09, 1e-10));
}

TEST_F(TwoDoubleVectorsForVectorUtils, CanBeFiltered) {  // NOLINT
  std::vector<double> filtered = util::vector_utils::GetEntriesAtIndices(double_vector_a_, {2});
  ASSERT_THAT(filtered.size(), 1);
  ASSERT_THAT(filtered[0], DoubleEq(8.9));
  ASSERT_THROW(util::vector_utils::GetEntriesAtIndices(double_vector_b_, {2, 3}), std::invalid_argument);
}
