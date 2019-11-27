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
                    transformation_matrix_({std::array<double, 4>({0.5, 0.0, sqrt(3) / 2.0, 0.0}),
                                            std::array<double, 4>({0.0, 1.0, 0.0, 0.0}),
                                            std::array<double, 4>({-sqrt(3) / 2.0, 0.0, 0.5, 0.0}),
                                            std::array<double, 4>({0.0, 0.0, 0.0, 1.0})}) {}

 protected:
  baf::ControlPoint control_point_a_;
  baf::ControlPoint control_point_b_;
  baf::ControlPoint control_point_c_;
  std::array<std::array<double, 4>, 4> transformation_matrix_;
};

bool operator==(baf::ControlPoint const &lhs, baf::ControlPoint const &rhs) {
  for (Dimension current_dimension{0}; current_dimension < Dimension{rhs.GetDimensionality()}; ++current_dimension) {
    if (!util::numeric_settings::AreEqual(rhs[current_dimension] , lhs[current_dimension])) return false;
  }
  return true;
}

TEST_F(AControlPoint, CanBeCopied) {  // NOLINT
  baf::ControlPoint copy_of_a(control_point_a_);
  ASSERT_TRUE(operator==(copy_of_a, control_point_a_));
}

TEST_F(AControlPoint, CanBeMoved) {  // NOLINT
  baf::ControlPoint copy_of_a(std::move(control_point_a_));
  ASSERT_TRUE(operator==(copy_of_a, control_point_a_));
}

TEST_F(AControlPoint, CanBeAssigned) {  // NOLINT
  baf::ControlPoint copy_of_a = control_point_a_;
  ASSERT_TRUE(operator==(copy_of_a, control_point_a_));
}

TEST_F(AControlPoint, CanBeMoveAssigned) {  // NOLINT
  baf::ControlPoint copy_of_a = std::move(control_point_a_);
  ASSERT_TRUE(operator==(copy_of_a, control_point_a_));
}

TEST_F(AControlPoint, ReturnsCorrectDimensionalityOf2) {  // NOLINT
  ASSERT_THAT(control_point_a_.GetDimensionality(), 2);
}

TEST_F(AControlPoint, Returns1ForDimension0) {  // NOLINT
  ASSERT_THAT(control_point_a_.GetValueForDimension(Dimension{0}), DoubleEq(1.0));
  ASSERT_THAT(control_point_a_[Dimension{0}], DoubleEq(1.0));
}

TEST_F(AControlPoint, Returns2ForDimension1) {  // NOLINT
  ASSERT_THAT(control_point_a_.GetValueForDimension(Dimension{1}), DoubleEq(2.0));
  ASSERT_THAT(control_point_a_[Dimension{1}], DoubleEq(2.0));
}

TEST_F(AControlPoint, Returns2_5ForDimension1AfterSetValue) {  // NOLINT
  baf::ControlPoint control_point(control_point_a_);
  control_point_a_.SetValue(Dimension{1}, 2.5);
  ASSERT_THAT(control_point_a_[Dimension{1}], DoubleEq(2.5));
}

TEST_F(AControlPoint, ReturnsCorrectResultAfterAdditionOfControlPoints) {  // NOLINT
  ASSERT_THAT((control_point_a_ + control_point_b_)[Dimension{0}], DoubleEq(3.0));
  ASSERT_THAT((control_point_a_ + control_point_b_)[Dimension{1}], DoubleEq(1.0));
}

TEST_F(AControlPoint, ReturnsCorrectResultAfterSubtractionOfControlPoints) {  // NOLINT
  ASSERT_THAT((control_point_a_ - control_point_b_)[Dimension{0}], DoubleEq(-1.0));
  ASSERT_THAT((control_point_a_ - control_point_b_)[Dimension{1}], DoubleEq(3.0));
}

TEST_F(AControlPoint, ReturnsCorrectResultAfterMultiplicationWithScalar) {  // NOLINT
  ASSERT_THAT((control_point_a_ * 0.5)[Dimension{0}], DoubleEq(0.5));
  ASSERT_THAT((control_point_a_ * 0.5)[Dimension{1}], DoubleEq(1.0));
  ASSERT_THAT((0.5 * control_point_a_)[Dimension{0}], DoubleEq(0.5));
  ASSERT_THAT((0.5 * control_point_a_)[Dimension{1}], DoubleEq(1.0));
}

TEST_F(AControlPoint, PerformsTransformationWithNoScalingCorrectly) {  // NOLINT
  baf::ControlPoint control_point_t = control_point_c_.Transform(transformation_matrix_, {1.0, 1.0, 1.0});
  ASSERT_THAT(control_point_t[Dimension{0}], DoubleEq((10 + 19 * sqrt(3)) / 20));
  ASSERT_THAT(control_point_t[Dimension{1}], DoubleEq(0.7));
  ASSERT_THAT(control_point_t[Dimension{2}], DoubleNear((19 - 10 * sqrt(3)) / 20, 1e-10));
}

TEST_F(AControlPoint, PerformsTransformationWithScalingOf0_5And1And0_7Correctly) {  // NOLINT
  baf::ControlPoint control_point_t = control_point_c_.Transform(transformation_matrix_, {0.5, 1.0, 0.7});
  ASSERT_THAT(control_point_t[Dimension{0}], DoubleEq((50 + 133 * sqrt(3)) / 200));
  ASSERT_THAT(control_point_t[Dimension{1}], DoubleEq(0.7));
  ASSERT_THAT(control_point_t[Dimension{2}], DoubleEq((133 - 50 * sqrt(3)) / 200));
}

TEST_F(AControlPoint, ReturnsEuclideanNormOfSquareRootOf5) {  // NOLINT
  ASSERT_THAT(control_point_a_.GetEuclideanNorm(), DoubleEq(sqrt(5.0)));
}
