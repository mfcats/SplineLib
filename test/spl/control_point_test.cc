/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include <functional>

#include "gmock/gmock.h"

#include "src/spl/control_point.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

using namespace splinelib::src;

class AControlPoint : public Test {
 public:
  AControlPoint() : control_point_a_({1.0, 2.0}), control_point_b_({2.0, -1.0}), control_point_c_({1.0, 0.7, 1.9}),
                    transformation_matrix_({std::array<double, 4>({0.5, 0.0, 0.866, 0.0}),
                                            std::array<double, 4>({0.0, 1.0, 0.0, 0.0}),
                                            std::array<double, 4>({-0.866, 0.0, 0.5, 0.0}),
                                            std::array<double, 4>({0.0, 0.0, 0.0, 1.0})}),
                    scaling_({1.0, 1.0, 1.0}) {}

 protected:
  baf::ControlPoint control_point_a_;
  baf::ControlPoint control_point_b_;
  baf::ControlPoint control_point_c_;
  std::array<std::array<double, 4>, 4> transformation_matrix_;
  std::array<double, 3> scaling_;
};

bool operator== (baf::ControlPoint const &lhs, baf::ControlPoint const &rhs) {
  for (Dimension i{0}; i < Dimension{rhs.GetDimensionality()}; ++i) {
    if (rhs[i] != lhs[i]) return false;
  }
  return true;
}

TEST_F(AControlPoint, CopyConstructorWorks) {  // NOLINT
  baf::ControlPoint copy_of_a(control_point_a_);
  ASSERT_TRUE(operator==(copy_of_a, control_point_a_));
}

TEST_F(AControlPoint, MoveConstructorWorks) {  // NOLINT
  baf::ControlPoint copy_of_a(std::move(control_point_a_));
  ASSERT_TRUE(operator==(copy_of_a, control_point_a_));
}

TEST_F(AControlPoint, CopyAssignmentOperatorWorks) {  // NOLINT
  baf::ControlPoint copy_of_a = control_point_a_;
  ASSERT_TRUE(operator==(copy_of_a, control_point_a_));
}

TEST_F(AControlPoint, MoveAssignmentOperatorWorks) {  // NOLINT
  baf::ControlPoint copy_of_a = std::move(control_point_a_);
  ASSERT_TRUE(operator==(copy_of_a, control_point_a_));
}

TEST_F(AControlPoint, ReturnsCorrectDimension) { // NOLINT
  ASSERT_THAT(control_point_a_.GetDimensionality(), 2);
}

TEST_F(AControlPoint, Returns1ForDimension0) { // NOLINT
  ASSERT_THAT(control_point_a_.GetValueForDimension(Dimension{0}), DoubleEq(1.0));
  ASSERT_THAT(control_point_a_[Dimension{0}], DoubleEq(1.0));
}

TEST_F(AControlPoint, Returns2ForDimension1) { // NOLINT
  ASSERT_THAT(control_point_a_.GetValueForDimension(Dimension{1}), DoubleEq(2.0));
  ASSERT_THAT(control_point_a_[Dimension{1}], DoubleEq(2.0));
}

TEST_F(AControlPoint, Returns2_5ForDimension1AfterSetValue) { // NOLINT
  baf::ControlPoint control_point(control_point_a_);
  control_point_a_.SetValue(1, 2.5);
  ASSERT_THAT(control_point_a_[Dimension{1}], DoubleEq(2.5));
}

TEST_F(AControlPoint, ReturnsCorrectResultAfterAdditionOfControlPoints) { // NOLINT
  ASSERT_THAT((control_point_a_ + control_point_b_)[Dimension{0}], DoubleEq(3.0));
  ASSERT_THAT((control_point_a_ + control_point_b_)[Dimension{1}], DoubleEq(1.0));
}

TEST_F(AControlPoint, ReturnsCorrectResultAfterSubtractionOfControlPoints) { // NOLINT
  ASSERT_THAT((control_point_a_ - control_point_b_)[Dimension{0}], DoubleEq(-1.0));
  ASSERT_THAT((control_point_a_ - control_point_b_)[Dimension{1}], DoubleEq(3.0));
}

TEST_F(AControlPoint, ReturnsCorrectResultAfterMultiplicationWithScalar) { // NOLINT
  ASSERT_THAT((control_point_a_ * 0.5)[Dimension{0}], DoubleEq(0.5));
  ASSERT_THAT((control_point_a_ * 0.5)[Dimension{1}], DoubleEq(1.0));
  ASSERT_THAT((0.5 * control_point_a_)[Dimension{0}], DoubleEq(0.5));
  ASSERT_THAT((0.5 * control_point_a_)[Dimension{1}], DoubleEq(1.0));
}

TEST_F(AControlPoint, PerformsTransformationCorrectly) { // NOLINT
  baf::ControlPoint control_point_t = control_point_c_.Transform(transformation_matrix_, scaling_);
  ASSERT_THAT(control_point_t[Dimension{0}], DoubleEq(2.1454));
}

TEST_F(AControlPoint, PerformsScalingTransformation) { // NOLINT
  std::array<double, 3> customScaling = {0.5, 1, 0.7};
  baf::ControlPoint control_point_t = control_point_c_.Transform(transformation_matrix_, customScaling);
  ASSERT_THAT(control_point_t[Dimension{0}], DoubleNear(1.40178, 0.00001));
}

TEST_F(AControlPoint, ReturnsEuclideanNormOf) { // NOLINT
  ASSERT_THAT(control_point_a_.GetEuclideanNorm(), DoubleEq(sqrt(5.0)));
}
