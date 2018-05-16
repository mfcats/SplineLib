/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <numeric_settings.h>
#include "gmock/gmock.h"

#include "nurbs.h"

using testing::Test;
using testing::DoubleNear;

class ARationalBSpline : public Test {
 public:
  ARationalBSpline() {
    std::array<baf::KnotVector, 1> knot_vector = {baf::KnotVector({0, 0, 0, 1, 2, 3, 3, 3})};
    std::array<int, 1> degree = {2};
    std::vector<double> weights = {1, 4, 1, 1, 1};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({5.0, -1.0}))
    };
    nurbs = std::make_unique<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::unique_ptr<spl::NURBS<1>> nurbs;
};

TEST_F(ARationalBSpline, ReturnsCorrectValues) {
  ASSERT_THAT(nurbs->Evaluate({1.0}, {0})[0], DoubleNear(7.0 / 5.0, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({1.0}, {1})[0], DoubleNear(6.0 / 5.0, util::NumericSettings<double>::kEpsilon()));
}



