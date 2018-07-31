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

#include "b_spline.h"
#include "projection.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class ABSpline2 : public Test {
 public:
  ABSpline2() {
    std::array<baf::KnotVector, 1> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.2}, ParamCoord{0.4},
                          ParamCoord{0.6}, ParamCoord{0.8}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1},
                          ParamCoord{1}})};
    std::array<Degree, 1> degree = {Degree{3}};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({100, 100})),
        baf::ControlPoint(std::vector<double>({140, 196})),
        baf::ControlPoint(std::vector<double>({200, 240})),
        baf::ControlPoint(std::vector<double>({260, 164})),
        baf::ControlPoint(std::vector<double>({340, 164})),
        baf::ControlPoint(std::vector<double>({400, 240})),
        baf::ControlPoint(std::vector<double>({460, 196})),
        baf::ControlPoint(std::vector<double>({500, 100}))
    };
    knot_vector_ptr[0] = std::make_shared<baf::KnotVector>(knot_vector[0]);
    b_spline = std::make_unique<spl::BSpline<1>>(knot_vector_ptr, degree, control_points);
  }

 protected:
  std::unique_ptr<spl::Spline<1>> b_spline;
  std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_ptr;
};

TEST_F(ABSpline2, CloseToCenter) { // NOLINT
  ASSERT_THAT(spl::Projection<1>::ProjectionOnSpline({332, 200}, std::move(b_spline))[0],
              DoubleNear(0.6223419238, 0.0001));
}

TEST_F(ABSpline2, RightOfLastKnot) { // NOLINT
  ASSERT_THAT(spl::Projection<1>::ProjectionOnSpline({800, 200}, std::move(b_spline))[0],
              DoubleEq(1));
}

TEST_F(ABSpline2, LeftOfFirstKnot) { // NOLINT
  ASSERT_THAT(spl::Projection<1>::ProjectionOnSpline({0, 0}, std::move(b_spline))[0],
              DoubleEq(0));
}
