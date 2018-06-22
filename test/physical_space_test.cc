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

#include "physical_space.h"
#include "weighted_physical_space.h"

using testing::Test;
using testing::DoubleEq;

class A1DPhysicalSpace : public Test {
 public:
  A1DPhysicalSpace() {
    control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({5.0, -1.0}))
    };
    physical_space = spl::PhysicalSpace<1>(control_points, {5});
  }

 protected:
  spl::PhysicalSpace<1> physical_space;
  std::vector<baf::ControlPoint> control_points;
};

TEST_F(A1DPhysicalSpace, ThrowsForDifferingGivenNumberOfControlPointsAndLengthOfControlPointVector) {
  ASSERT_THROW(spl::PhysicalSpace<1>(control_points, {4}), std::runtime_error);
}

TEST_F(A1DPhysicalSpace, ThrowsForDifferingDimensionsOfControlPoints) {
  control_points.emplace_back(std::vector<double>({0.0}));
  ASSERT_THROW(spl::PhysicalSpace<1>(control_points, {6}), std::runtime_error);
}

TEST_F(A1DPhysicalSpace, ReturnsCorrectFirstControlPoint) {
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{0}).GetValue(0), DoubleEq(0.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{0}).GetValue(1), DoubleEq(0.0));
}

TEST_F(A1DPhysicalSpace, ReturnsCorrectInnerControlPoint) {
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{2}).GetValue(0), DoubleEq(3.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{2}).GetValue(1), DoubleEq(2.0));
}

TEST_F(A1DPhysicalSpace, ReturnsCorrectLastControlPoint) {
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{4}).GetValue(0), DoubleEq(5.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{4}).GetValue(1), DoubleEq(-1.0));
}

class A1DWeightedPhysicalSpace : public A1DPhysicalSpace {
 public:
  A1DWeightedPhysicalSpace() {
    weights_ = {0.5, 0.75, 0.8, 1.0, 1.2};
    weighted_physical_space = spl::WeightedPhysicalSpace<1>(control_points, weights_, {5});
  }

 protected:
  spl::WeightedPhysicalSpace<1> weighted_physical_space;
  std::vector<double> weights_;
};

TEST_F(A1DWeightedPhysicalSpace, ThrowsForDifferingNumberOfControlPointsAndWeights) {
  control_points.emplace_back(std::vector<double>({0.0}));
  ASSERT_THROW(spl::WeightedPhysicalSpace<1>(control_points, weights_, {6}), std::runtime_error);
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectWeight) {
  ASSERT_THAT(weighted_physical_space.GetWeight({2}), DoubleEq(0.8));
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectFirstHomogenousControlPoint) {
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{0}).GetValue(0), DoubleEq(0.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{0}).GetValue(1), DoubleEq(0.0));
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectInnerHomogenousControlPoint) {
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{2}).GetValue(0), DoubleEq(2.4));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{2}).GetValue(1), DoubleEq(1.6));
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectLastHomogenousControlPoint) {
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{4}).GetValue(0), DoubleEq(6.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(std::array<int, 1>{4}).GetValue(1), DoubleEq(-1.2));
}

class A2DPhysicalSpace : public Test {
 public:
  A2DPhysicalSpace() {
    control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.5, 2.5})),
        baf::ControlPoint(std::vector<double>({5.0, 1.0}))
    };
    physical_space = spl::PhysicalSpace<2>(control_points, {3, 2});
  }

 protected:
  spl::PhysicalSpace<2> physical_space;
  std::vector<baf::ControlPoint> control_points;
};

TEST_F(A2DPhysicalSpace, ThrowsForDifferingGivenNumberOfControlPointsAndLengthOfControlPointVector) {
  ASSERT_THROW(spl::PhysicalSpace<2>(control_points, {3, 1}), std::runtime_error);
}

TEST_F(A2DPhysicalSpace, ThrowsForDifferingDimensionsOfControlPoints) {
  control_points.emplace_back(std::vector<double>({0.0}));
  ASSERT_THROW(spl::PhysicalSpace<2>(control_points, {1, 7}), std::runtime_error);
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectFirstControlPoint) {
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{0, 0}).GetValue(0), DoubleEq(0.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{0, 0}).GetValue(1), DoubleEq(0.0));
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectInnerControlPoint) {
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{1, 1}).GetValue(0), DoubleEq(1.5));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{1, 1}).GetValue(1), DoubleEq(2.5));
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectLastControlPoint) {
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{2, 1}).GetValue(0), DoubleEq(5.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{2, 1}).GetValue(1), DoubleEq(1.0));
}
