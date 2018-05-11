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

using testing::Test;
using testing::DoubleEq;

class ARationalBSpline : public Test {
 public:
  ARationalBSpline() {
    std::array<KnotVector, 1> knot_vector = {KnotVector({0, 0, 0, 1, 2, 3, 3, 3})};
    std::array<int, 1> degree = {2};
    std::vector<double> weights = {1, 4, 1, 1, 1};
    std::vector<ControlPoint> control_points = {
        ControlPoint(std::vector<double>({0.0, 0.0})),
        ControlPoint(std::vector<double>({1.0, 1.0})),
        ControlPoint(std::vector<double>({3.0, 2.0})),
        ControlPoint(std::vector<double>({4.0, 1.0})),
        ControlPoint(std::vector<double>({5.0, -1.0}))
    };
    nurbs = std::make_unique<NURBS < 1>>(knot_vector, degree, weights, control_points);
  }

 protected:
  std::unique_ptr<NURBS < 1>> nurbs;
};

TEST_F(ARationalBSpline,) {
  ASSERT_NEAR(nurbs->Evaluate({1.0}, {0})[0], 3.5, 0.00005);
  ASSERT_NEAR(nurbs->Evaluate({1.0}, {1})[0], 3.0, 0.00005);
  ASSERT_NEAR(nurbs->Evaluate({1.0}, {2})[0], 2.5, 0.00005);
}



