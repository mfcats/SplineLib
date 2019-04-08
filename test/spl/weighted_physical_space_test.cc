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

#include "weighted_physical_space.h"

using testing::Test;
using testing::DoubleEq;

class A1DWeightedPhysicalSpace : public Test {
 public:
  A1DWeightedPhysicalSpace() {
    control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({5.0, -1.0}))
    };
    weights_ = {0.5, 0.75, 0.8, 1.0, 1.2};
    weighted_physical_space = spl::WeightedPhysicalSpace<1>(control_points, weights_, {5});
  }

 protected:
  spl::WeightedPhysicalSpace<1> weighted_physical_space;
  std::vector<baf::ControlPoint> control_points;
  std::vector<double> weights_;
};

TEST_F(A1DWeightedPhysicalSpace, ThrowsForDifferingNumberOfControlPointsAndWeights) {  // NOLINT
  control_points.emplace_back(std::vector<double>({0.0, 0.0}));
  ASSERT_THROW(spl::WeightedPhysicalSpace<1>(control_points, weights_, {6}), std::runtime_error);
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectWeight) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight({2}), DoubleEq(0.8));
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectFirstHomogenousControlPoint) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{0}).GetValue(0), DoubleEq(0.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{0}).GetValue(1), DoubleEq(0.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{0}).GetValue(2), DoubleEq(0.5));
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectInnerHomogenousControlPoint) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{2}).GetValue(0), DoubleEq(2.4));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{2}).GetValue(1), DoubleEq(1.6));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{2}).GetValue(2), DoubleEq(0.8));
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectLastHomogenousControlPoint) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{4}).GetValue(0), DoubleEq(6.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{4}).GetValue(1), DoubleEq(-1.2));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{4}).GetValue(2), DoubleEq(1.2));
}

class A2DWeightedPhysicalSpace : public Test {
 public:
  A2DWeightedPhysicalSpace() {
    control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.5, 2.5})),
        baf::ControlPoint(std::vector<double>({5.0, 1.0}))
    };
    weights_ = {0.5, 0.75, 0.8, 1.0, 1.2, 3.8};
    weighted_physical_space = spl::WeightedPhysicalSpace<2>(control_points, weights_, {3, 2});
  }

 protected:
  spl::WeightedPhysicalSpace<2> weighted_physical_space;
  std::vector<baf::ControlPoint> control_points;
  std::vector<double> weights_;
};

TEST_F(A2DWeightedPhysicalSpace, ThrowsForDifferingNumberOfControlPointsAndWeights) {  // NOLINT
  weights_.emplace_back(1.0);
  ASSERT_THROW(spl::WeightedPhysicalSpace<2>(control_points, weights_, {3, 2}), std::runtime_error);
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectFirstWeightFor2DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight({0, 0}), DoubleEq(0.5));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectInnerWeightFor2DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight({0, 1}), DoubleEq(1.0));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectLastWeightFor2DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight({2, 1}), DoubleEq(3.8));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectFirstWeightFor1DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight({0}), DoubleEq(0.5));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectInnerWeightFor1DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight({3}), DoubleEq(1.0));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectLastWeightFor1DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight({5}), DoubleEq(3.8));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectFirstHomogenousControlPointFor2DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{0, 0}).GetValue(0), DoubleEq(0.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{0, 0}).GetValue(1), DoubleEq(0.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{0, 0}).GetValue(2), DoubleEq(0.5));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectInnerHomogenousControlPointFor2DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{1, 1}).GetValue(0), DoubleEq(1.8));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{1, 1}).GetValue(1), DoubleEq(3.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{1, 1}).GetValue(2), DoubleEq(1.2));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectLastHomogenousControlPointFor2DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{2, 1}).GetValue(0), DoubleEq(19.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{2, 1}).GetValue(1), DoubleEq(3.8));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{2, 1}).GetValue(2), DoubleEq(3.8));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectFirstHomogenousControlPointFor1DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{0}).GetValue(0), DoubleEq(0.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{0}).GetValue(1), DoubleEq(0.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{0}).GetValue(2), DoubleEq(0.5));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectInnerHomogenousControlPointFor1DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{4}).GetValue(0), DoubleEq(1.8));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{4}).GetValue(1), DoubleEq(3.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{4}).GetValue(2), DoubleEq(1.2));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectLastHomogenousControlPointFor1DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{5}).GetValue(0), DoubleEq(19.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{5}).GetValue(1), DoubleEq(3.8));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 2>{5}).GetValue(2), DoubleEq(3.8));
}

/* 1-dimensional nurbs spline with following properties :
 * KnotVector = {0, 0, 0, 1, 2, 3, 3, 3}
 * ControlPoints = {{0, 0}, {1, 1}, {3, 2}, {4, 1}, {5, -1}}
 * Weights = {1, 4, 1, 1, 1}
*/

class A1DWeightedPhysicalSpace_a : public Test {
 public:
  A1DWeightedPhysicalSpace_a() {
    control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({5.0, -1.0}))
    };
    weights_ = {1, 4, 1, 1, 1};
    weighted_physical_space = spl::WeightedPhysicalSpace<1>(control_points, weights_, {5});
  }

 protected:
  spl::WeightedPhysicalSpace<1> weighted_physical_space;
  std::vector<baf::ControlPoint> control_points;
  std::vector<double> weights_;
};

TEST_F(A1DWeightedPhysicalSpace_a, ReturnsHomogeneousControlPoint) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{1}).GetValue(0), DoubleEq(4.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{1}).GetValue(1), DoubleEq(4.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{1}).GetValue(2), DoubleEq(4.0));
}
