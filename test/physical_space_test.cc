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

TEST_F(A1DPhysicalSpace, ThrowsForDifferingGivenNumberOfControlPointsAndLengthOfControlPointVector) {  // NOLINT
  ASSERT_THROW(spl::PhysicalSpace<1>(control_points, {4}), std::runtime_error);
}

TEST_F(A1DPhysicalSpace, ThrowsForDifferingDimensionsOfControlPoints) {  // NOLINT
  control_points.emplace_back(std::vector<double>({0.0}));
  ASSERT_THROW(spl::PhysicalSpace<1>(control_points, {6}), std::runtime_error);
}

TEST_F(A1DPhysicalSpace, ReturnsCorrectFirstControlPoint) {  // NOLINT
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{0}).GetValue(0), DoubleEq(0.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{0}).GetValue(1), DoubleEq(0.0));
}

TEST_F(A1DPhysicalSpace, ReturnsCorrectInnerControlPoint) {  // NOLINT
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{2}).GetValue(0), DoubleEq(3.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{2}).GetValue(1), DoubleEq(2.0));
}

TEST_F(A1DPhysicalSpace, ReturnsCorrectLastControlPoint) {  // NOLINT
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{4}).GetValue(0), DoubleEq(5.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{4}).GetValue(1), DoubleEq(-1.0));
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

TEST_F(A2DPhysicalSpace, ThrowsForDifferingGivenNumberOfControlPointsAndLengthOfControlPointVector) {  // NOLINT
  ASSERT_THROW(spl::PhysicalSpace<2>(control_points, {3, 1}), std::runtime_error);
}

TEST_F(A2DPhysicalSpace, ThrowsForDifferingDimensionsOfControlPoints) {  // NOLINT
  control_points.emplace_back(std::vector<double>({0.0}));
  ASSERT_THROW(spl::PhysicalSpace<2>(control_points, {1, 7}), std::runtime_error);
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectFirstControlPointFor2DIndex) {  // NOLINT
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{0, 0}).GetValue(0), DoubleEq(0.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{0, 0}).GetValue(1), DoubleEq(0.0));
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectInnerControlPointFor2DIndex) {  // NOLINT
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{1, 1}).GetValue(0), DoubleEq(1.5));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{1, 1}).GetValue(1), DoubleEq(2.5));
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectLastControlPointFor2DIndex) {  // NOLINT
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{2, 1}).GetValue(0), DoubleEq(5.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{2, 1}).GetValue(1), DoubleEq(1.0));
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectFirstControlPointFor1DIndex) {  // NOLINT
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{0}).GetValue(0), DoubleEq(0.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{0}).GetValue(1), DoubleEq(0.0));
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectInnerControlPointFor1DIndex) {  // NOLINT
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{4}).GetValue(0), DoubleEq(1.5));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{4}).GetValue(1), DoubleEq(2.5));
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectLastControlPointFor1DIndex) {  // NOLINT
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{5}).GetValue(0), DoubleEq(5.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{5}).GetValue(1), DoubleEq(1.0));
}
