/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "b_spline.h"

#include "gmock/gmock.h"

using testing::Test;
using testing::DoubleEq;

class A2DBSpline : public Test {
 public:
  A2DBSpline() {
    std::array<KnotVector, 2> knot_vector = {KnotVector({0, 0, 0, 1, 1, 1}), KnotVector({0, 0, 0, 1, 1, 1})};
    std::array<int, 2> degree = {2, 2};
    std::vector<ControlPoint> control_points = {
        ControlPoint(std::vector<double>({-1.0, -1.0, 0.0})),
        ControlPoint(std::vector<double>({0.0, -1.0, 0.0})),
        ControlPoint(std::vector<double>({1.0, -1.0, 0.0})),
        ControlPoint(std::vector<double>({-1.0, 0.0, 0.0})),
        ControlPoint(std::vector<double>({0.0, 0.0, 1.0})),
        ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
        ControlPoint(std::vector<double>({-1.0, 1.0, 0.0})),
        ControlPoint(std::vector<double>({0.0, 1.0, 0.0})),
        ControlPoint(std::vector<double>({1.0, 1.0, 0.0}))
    };
    b_spline = std::make_unique<BSpline<2>>(knot_vector, degree, control_points);
  }

 protected:
  std::unique_ptr<BSpline<2>> b_spline;
};

TEST_F(A2DBSpline, Corner) {
  ASSERT_NEAR(b_spline->Evaluate({0.0, 0.0}, {0})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({0.0, 0.0}, {1})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({0.0, 0.0}, {2})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, EdgeDim0) {
  ASSERT_NEAR(b_spline->Evaluate({0.0, 0.33333}, {0})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({0.0, 0.33333}, {1})[0], -0.33333, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({0.0, 0.33333}, {2})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, EdgeDim1) {
  ASSERT_NEAR(b_spline->Evaluate({0.33333, 0.0}, {0})[0], -0.33333, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({0.33333, 0.0}, {1})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({0.33333, 0.0}, {2})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, Center) {
  ASSERT_NEAR(b_spline->Evaluate({0.5, 0.5}, {0})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({0.5, 0.5}, {1})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({0.5, 0.5}, {2})[0], 0.25, 0.00005);
}

TEST_F(A2DBSpline, Random) {
  ASSERT_NEAR(b_spline->Evaluate({0.75, 0.25}, {0})[0], 0.5, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({0.75, 0.25}, {1})[0], -0.5, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({0.75, 0.25}, {2})[0], 0.14063, 0.00005);
}