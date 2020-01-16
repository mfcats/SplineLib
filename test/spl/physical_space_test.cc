/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/spl/physical_space.h"

using testing::Test;
using testing::DoubleEq;

using namespace splinelib::src;

class A1DPhysicalSpace : public Test {
 public:
  A1DPhysicalSpace() {
    control_points = {
        spl::ControlPoint(std::vector<double>({0.0, 0.0})),
        spl::ControlPoint(std::vector<double>({1.0, 1.0})),
        spl::ControlPoint(std::vector<double>({3.0, 2.0})),
        spl::ControlPoint(std::vector<double>({4.0, 1.0})),
        spl::ControlPoint(std::vector<double>({5.0, -1.0}))
    };
    physical_space = spl::PhysicalSpace<1>(control_points, {5});
  }

 protected:
  spl::PhysicalSpace<1> physical_space;
  std::vector<spl::ControlPoint> control_points;
};

TEST_F(A1DPhysicalSpace, ThrowsForDifferingGivenNumberOfControlPointsAndLengthOfControlPointVector) {  // NOLINT
  ASSERT_THROW(spl::PhysicalSpace<1>(control_points, {4}), std::runtime_error);
}

TEST_F(A1DPhysicalSpace, ThrowsForDifferingDimensionsOfControlPoints) {  // NOLINT
  control_points.emplace_back(std::vector<double>({0.0}));
  ASSERT_THROW(spl::PhysicalSpace<1>(control_points, {6}), std::runtime_error);
}

TEST_F(A1DPhysicalSpace, ReturnsCorrectControlPoint) {  // NOLINT
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{2})[Dimension{0}], DoubleEq(3.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 1>{2})[Dimension{1}], DoubleEq(2.0));
}

TEST_F(A1DPhysicalSpace, ReturnsCorrectNumberOfPointsInEachDirection) {  // NOLINT
  ASSERT_THAT(physical_space.GetNumberOfPointsPerDirection()[0], 5);
}

TEST_F(A1DPhysicalSpace, ReturnsCorrectMaximumPointIndexForEachDirection) {  // NOLINT
  ASSERT_THAT(physical_space.GetMaximumPointIndexPerDirection()[0], 4);
}

TEST_F(A1DPhysicalSpace, ReturnsCorrectDimension) {  // NOLINT
  ASSERT_THAT(physical_space.GetDimensionality(), 2);
}

TEST_F(A1DPhysicalSpace, AddsAndSetsNewControlPoint) {  // NOLINT
  ASSERT_THAT(physical_space.GetTotalNumberOfControlPoints(), 5);
  // TODO(Corinna, Christoph): It should not be possible to add a control point without specifying a dimension.
  // physical_space.AddControlPoints(1);
  // physical_space.SetControlPoint(5, spl::ControlPoint(std::vector<double>({6.0, 0.0})));
  // ASSERT_THAT(physical_space.GetTotalNumberOfControlPoints(), 6);
  // ASSERT_THAT(physical_space.GetControlPoint(5)[0], DoubleEq(6.0));
  // ASSERT_THAT(physical_space.GetControlPoint(5)[1], DoubleEq(0.0));
  physical_space.SetControlPoint(4, spl::ControlPoint(std::vector<double>({6.0, 0.0})));
  ASSERT_THAT(physical_space.GetControlPoint(4)[Dimension{0}], DoubleEq(6.0));
  ASSERT_THAT(physical_space.GetControlPoint(4)[Dimension{1}], DoubleEq(0.0));
}

TEST_F(A1DPhysicalSpace, RemovesControlPoint) {  // NOLINT
  ASSERT_THAT(physical_space.GetTotalNumberOfControlPoints(), 5);
  // TODO(Corinna, Christoph): Same problem as above.
  physical_space.RemoveControlPoints(2);
  ASSERT_THAT(physical_space.GetTotalNumberOfControlPoints(), 3);
  ASSERT_THAT(physical_space.GetControlPoint(2)[Dimension{0}], DoubleEq(3.0));
  ASSERT_THAT(physical_space.GetControlPoint(2)[Dimension{1}], DoubleEq(2.0));
}

TEST_F(A1DPhysicalSpace, ReturnsCorrectExpansion) {  // NOLINT
  ASSERT_THAT(physical_space.GetExpansion(), DoubleEq(5.0));
}

TEST_F(A1DPhysicalSpace, ReturnsCorrectMaximumDistanceFromOrigin) {  // NOLINT
  ASSERT_THAT(physical_space.GetMaximumDistanceFromOrigin(), DoubleEq(sqrt(26)));
}

TEST_F(A1DPhysicalSpace, ReturnsCorrectDividedControlPoints) {  // NOLINT
  ASSERT_THAT(physical_space.SplitControlPoints(Dimension{0}, 2, 2).size(), 2);
  ASSERT_THAT(physical_space.SplitControlPoints(Dimension{0}, 2, 2)[0][Dimension{0}], DoubleEq(3.0));
}

TEST_F(A1DPhysicalSpace, EvaluatesThatCopiedSpaceEqualsOriginalSpace) {  // NOLINT
  spl::PhysicalSpace<1> copy(physical_space);
  ASSERT_THAT(physical_space.AreEqual(copy), true);
}

class A2DPhysicalSpace : public Test {
 public:
  A2DPhysicalSpace() {
    control_points = {
        spl::ControlPoint(std::vector<double>({0.0, 0.0})),
        spl::ControlPoint(std::vector<double>({1.0, 1.0})),
        spl::ControlPoint(std::vector<double>({3.0, 2.0})),
        spl::ControlPoint(std::vector<double>({0.0, 2.0})),
        spl::ControlPoint(std::vector<double>({1.5, 2.5})),
        spl::ControlPoint(std::vector<double>({5.0, 1.0}))
    };
    physical_space = spl::PhysicalSpace<2>(control_points, {3, 2});
  }

 protected:
  spl::PhysicalSpace<2> physical_space;
  std::vector<spl::ControlPoint> control_points;
};

TEST_F(A2DPhysicalSpace, ThrowsForDifferingGivenNumberOfControlPointsAndLengthOfControlPointVector) {  // NOLINT
  ASSERT_THROW(spl::PhysicalSpace<2>(control_points, {3, 1}), std::runtime_error);
}

TEST_F(A2DPhysicalSpace, ThrowsForDifferingDimensionsOfControlPoints) {  // NOLINT
  control_points.emplace_back(std::vector<double>({0.0}));
  ASSERT_THROW(spl::PhysicalSpace<2>(control_points, {1, 7}), std::runtime_error);
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectLastControlPointFor2DIndex) {  // NOLINT
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{2, 1})[Dimension{0}],
              DoubleEq(5.0));
  ASSERT_THAT(physical_space.GetControlPoint(std::array<int, 2>{2, 1})[Dimension{1}],
              DoubleEq(1.0));
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectInnerControlPointFor1DIndex) {  // NOLINT
  ASSERT_THAT(physical_space.GetControlPoint(4)[Dimension{0}], DoubleEq(1.5));
  ASSERT_THAT(physical_space.GetControlPoint(4)[Dimension{1}], DoubleEq(2.5));
}

TEST_F(A2DPhysicalSpace, ReturnsDefaultWeight) {  // NOLINT
  ASSERT_THAT(physical_space.GetWeight(std::array<int, 2>{1, 2}).Get(), DoubleEq(1.0));
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectNumberOfPointsInEachDirection) {  // NOLINT
  ASSERT_THAT(physical_space.GetNumberOfPointsPerDirection()[0], 3);
  ASSERT_THAT(physical_space.GetNumberOfPointsPerDirection()[1], 2);
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectMaximumPointIndexForEachDirection) {  // NOLINT
  ASSERT_THAT(physical_space.GetMaximumPointIndexPerDirection()[0], 2);
  ASSERT_THAT(physical_space.GetMaximumPointIndexPerDirection()[1], 1);
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectDimension) {  // NOLINT
  ASSERT_THAT(physical_space.GetDimensionality(), 2);
}

TEST_F(A2DPhysicalSpace, AddsAndSetsNewControlPoint) {  // NOLINT
  ASSERT_THAT(physical_space.GetTotalNumberOfControlPoints(), 6);
  // TODO(Corinna, Christoph): It should not be possible to add a control point without specifying a dimension.
  // physical_space.AddControlPoints(1);
  // physical_space.SetControlPoint(6, spl::ControlPoint(std::vector<double>({6.0, 0.0})));
  // ASSERT_THAT(physical_space.GetTotalNumberOfControlPoints(), 7);
  // ASSERT_THAT(physical_space.GetControlPoint(6)[0], DoubleEq(6.0));
  // ASSERT_THAT(physical_space.GetControlPoint(6)[1], DoubleEq(0.0));
  physical_space.SetControlPoint(5, spl::ControlPoint(std::vector<double>({6.0, 0.0})));
  ASSERT_THAT(physical_space.GetControlPoint(5)[Dimension{0}], DoubleEq(6.0));
  ASSERT_THAT(physical_space.GetControlPoint(5)[Dimension{1}], DoubleEq(0.0));
}

TEST_F(A2DPhysicalSpace, RemovesControlPoint) {  // NOLINT
  ASSERT_THAT(physical_space.GetTotalNumberOfControlPoints(), 6);
  // TODO(Corinna, Christoph): Same problem as above.
  // physical_space.RemoveControlPoints(2);
  // ASSERT_THAT(physical_space.GetTotalNumberOfControlPoints(), 4);
  ASSERT_THAT(physical_space.GetControlPoint(3)[Dimension{0}], DoubleEq(0.0));
  ASSERT_THAT(physical_space.GetControlPoint(3)[Dimension{1}], DoubleEq(2.0));
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectExpansion) {  // NOLINT
  ASSERT_THAT(physical_space.GetExpansion(), DoubleEq(5.0));
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectMaximumDistanceFromOrigin) {  // NOLINT
  ASSERT_THAT(physical_space.GetMaximumDistanceFromOrigin(), DoubleEq(sqrt(26)));
}

TEST_F(A2DPhysicalSpace, ReturnsCorrectDividedControlPoints) {  // NOLINT
  ASSERT_THAT(physical_space.SplitControlPoints(Dimension{0}, 1, 2).size(), 4);
  ASSERT_THAT(physical_space.SplitControlPoints(Dimension{0}, 1, 2)[2][Dimension{0}], DoubleEq(1.5));
}

TEST_F(A2DPhysicalSpace, EvaluatesThatCopiedSpaceEqualsOriginalSpace) {  // NOLINT
  spl::PhysicalSpace<2> copy(physical_space);
  ASSERT_THAT(physical_space.AreEqual(copy), true);
}
