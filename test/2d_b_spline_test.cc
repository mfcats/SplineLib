/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "b_spline_2_d.h"

#include "gmock/gmock.h"

using testing::Test;
using testing::DoubleEq;

class A2DBSpline : public Test {
 public:

    A2DBSpline() {
      std::vector<KnotVector> knot_tensor = {KnotVector({0, 0, 0, 0.5, 1, 1, 1}), KnotVector({0, 0, 0, 0.5, 1, 1, 1})};
      std::vector<int> degrees = {2, 2};
      std::vector<std::vector<ControlPoint>> control_points;

      std::vector<ControlPoint> control_points1 = {
          ControlPoint(std::vector<double>({-1.0, 1.0})),
          ControlPoint(std::vector<double>({0.0, -1.0})),
          ControlPoint(std::vector<double>({1.0, -1.0})),
      };
      std::vector<ControlPoint> control_points2 = {
          ControlPoint(std::vector<double>({-1.0, 0.0})),
          ControlPoint(std::vector<double>({0.0, 0.0})),
          ControlPoint(std::vector<double>({1.0, 0.0})),
      };
      std::vector<ControlPoint> control_points3 = {
          ControlPoint(std::vector<double>({-1.0, 1.0})),
          ControlPoint(std::vector<double>({0.0, 1.0})),
          ControlPoint(std::vector<double>({1.0, 1.0})),
      };

      control_points.push_back(control_points1);
      control_points.push_back(control_points2);
      control_points.push_back(control_points3);

      b_spline = std::make_unique<BSpline2D>(knot_tensor, degrees, control_points);
    }

 protected:
  std::unique_ptr<BSpline2D> b_spline;
};

TEST_F(A2DBSpline, Returns0_0For0AndDim0) {
  ASSERT_THAT(b_spline->Evaluate({0.0, 0.0}, {0})[0], DoubleEq(-1.0));
}

TEST_F(A2DBSpline, Returns0_0For0AndDim1) {
  ASSERT_THAT(b_spline->Evaluate({0.0, 0.0}, {1})[0], DoubleEq(-1.0));
}


/*
TEST_F(A2DBSpline, Returns4_0For5AndDim0) {
  ASSERT_THAT(b_spline->Evaluate(5.0, {0})[0], DoubleEq(4.0));
}

TEST_F(A2DBSpline, Returns0_0For5AndDim1) {
  ASSERT_THAT(b_spline->Evaluate(5.0, {1})[0], DoubleEq(0.0));
}

TEST_F(A2DBSpline, Returns1_5For2_5AndDim0) {
  ASSERT_THAT(b_spline->Evaluate(2.5, {0})[0], DoubleEq(1.5));
}

TEST_F(A2DBSpline, Returns0_0For0_0Dim0AndDer1) {
  ASSERT_THAT(b_spline->EvaluateDerivative(0.0, {0}, 1)[0], DoubleEq(0.0));
}

TEST_F(A2DBSpline, Returns1_0For0_0Dim1AndDer1) {
  ASSERT_THAT(b_spline->EvaluateDerivative(0.0, {1}, 1)[0], DoubleEq(2.0));
}

TEST_F(A2DBSpline, Returns12_0For5_0Dim0AndDer1) {
  ASSERT_THAT(b_spline->EvaluateDerivative(5.0, {0}, 1)[0], DoubleEq(0.0));
}

TEST_F(A2DBSpline, Returns0_325For2_25Dim1AndDer1) {
  ASSERT_THAT(b_spline->EvaluateDerivative(2.25, {1}, 1)[0], DoubleEq(0.325));
}
 */
