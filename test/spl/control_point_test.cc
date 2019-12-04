/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/spl/control_point.h"

using testing::DoubleEq;
using testing::Eq;
using testing::Test;

using namespace splinelib::src;

class ThreeControlPointsOfDimensionality2And3 : public Test {
 public:
  ThreeControlPointsOfDimensionality2And3()
    : control_point_a_({1.0, 0.7}), control_point_b_({2.0, -1.0}), control_point_c_({1.0, 0.7, 1.9}) {}

 protected:
  spl::ControlPoint control_point_a_;
  spl::ControlPoint control_point_b_;
  spl::ControlPoint control_point_c_;
};

TEST_F(ThreeControlPointsOfDimensionality2And3, ControlPointACanBeCopied) {  // NOLINT
  spl::ControlPoint copied_control_point_a(control_point_a_);
  ASSERT_THAT(copied_control_point_a, Eq(control_point_a_));
}

TEST_F(ThreeControlPointsOfDimensionality2And3, ControlPointACanBeMoved) {  // NOLINT
  spl::ControlPoint moved_control_point_a(std::move(control_point_a_));
  ASSERT_THAT(moved_control_point_a.GetValueForDimension(Dimension{1}), DoubleEq(0.7));
}

TEST_F(ThreeControlPointsOfDimensionality2And3, ControlPointACanBeAssigned) {  // NOLINT
  spl::ControlPoint control_point_to_assign(3);
  control_point_to_assign = control_point_a_;
  ASSERT_THAT(control_point_to_assign, Eq(control_point_a_));
}

TEST_F(ThreeControlPointsOfDimensionality2And3, ControlPointACanBeMoveAssigned) {  // NOLINT
  spl::ControlPoint control_point_to_move_assign(3);
  control_point_to_move_assign = std::move(control_point_a_);
  ASSERT_THAT(control_point_to_move_assign.GetValueForDimension(Dimension{1}), DoubleEq(0.7));
}

TEST_F(ThreeControlPointsOfDimensionality2And3, ControlPointAReturnsCorrectDimensionalityOf2) {  // NOLINT
  ASSERT_THAT(control_point_a_.GetDimensionality(), 2);
}

TEST_F(ThreeControlPointsOfDimensionality2And3, ControlPointAReturns1ForDimension0) {  // NOLINT
  ASSERT_THAT(control_point_a_.GetValueForDimension(Dimension{0}), DoubleEq(1.0));
  ASSERT_THAT(control_point_a_[Dimension{0}], DoubleEq(1.0));
}

TEST_F(ThreeControlPointsOfDimensionality2And3, ControlPointAReturns0_7ForDimension1) {  // NOLINT
  ASSERT_THAT(control_point_a_.GetValueForDimension(Dimension{1}), DoubleEq(0.7));
  ASSERT_THAT(control_point_a_[Dimension{1}], DoubleEq(0.7));
}

TEST_F(ThreeControlPointsOfDimensionality2And3, ControlPointAReturns2_5ForDimension1AfterSetValue) {  // NOLINT
  spl::ControlPoint control_point(control_point_a_);
  control_point_a_.SetValue(Dimension{1}, 2.5);
  ASSERT_THAT(control_point_a_[Dimension{1}], DoubleEq(2.5));
}

TEST_F(ThreeControlPointsOfDimensionality2And3, Return3AndMinus0_3AfterAdditionOfControlPointsAAndB) {  // NOLINT
  ASSERT_THAT((control_point_a_ + control_point_b_), Eq(spl::ControlPoint({3.0, -0.3})));
}

TEST_F(ThreeControlPointsOfDimensionality2And3, ThrowForAdditionOfControlPointsWithDimensionality2And3) {  //  NOLINT
  ASSERT_THROW((control_point_a_ + control_point_c_), std::invalid_argument);
}

TEST_F(ThreeControlPointsOfDimensionality2And3, ReturnMinus1And1_7AfterSubtractionOfControlPointsAAndB) {  // NOLINT
  ASSERT_THAT((control_point_a_ - control_point_b_), Eq(spl::ControlPoint({-1.0, 1.7})));
}

TEST_F(ThreeControlPointsOfDimensionality2And3, ThrowForSubtractionOfControlPointsAAndCWithDimensionality2And3) {  //  NOLINT
  ASSERT_THROW((control_point_a_ - control_point_c_), std::invalid_argument);
}

TEST_F(ThreeControlPointsOfDimensionality2And3, Return0_5And0_35AfterMultiplicationOfControlPointAWithScalar0_5) {  // NOLINT
  ASSERT_THAT((control_point_a_ * 0.5), Eq(spl::ControlPoint({0.5, 0.35})));
  ASSERT_THAT((0.5 * control_point_a_), Eq(spl::ControlPoint({0.5, 0.35})));
}

TEST_F(ThreeControlPointsOfDimensionality2And3, ControlPointBReturnsTwoNormOfSquareRootOf5) {  // NOLINT
  ASSERT_THAT(control_point_b_.ComputeTwoNorm(), DoubleEq(sqrt(5.0)));
}
