/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <any>

#include "gmock/gmock.h"

#include "b_spline.h"
#include "projection.h"
#include "writer.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class ABSplineSurface : public Test {
 public:
  ABSplineSurface() {
    std::array<baf::KnotVector, 2> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1},
                          ParamCoord{1}, ParamCoord{1}}),
         baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.25}, ParamCoord{0.5},
                          ParamCoord{0.75}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})};
    std::array<Degree, 2> degree = {Degree{3}, Degree{3}};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({-236, -197, -22})),
        baf::ControlPoint(std::vector<double>({-206, -117, -22})),
        baf::ControlPoint(std::vector<double>({-216, -27, 8})),
        baf::ControlPoint(std::vector<double>({-246, 62, -22})),
        baf::ControlPoint(std::vector<double>({-156, -177, 8})),
        baf::ControlPoint(std::vector<double>({-176, -97, 38})),
        baf::ControlPoint(std::vector<double>({-157, 20, 126})),
        baf::ControlPoint(std::vector<double>({-186, 142, 8})),
        baf::ControlPoint(std::vector<double>({-86, -157, 8})),
        baf::ControlPoint(std::vector<double>({-138, -113, -146})),
        baf::ControlPoint(std::vector<double>({-104, 14, -60})),
        baf::ControlPoint(std::vector<double>({-96, 102, 8})),
        baf::ControlPoint(std::vector<double>({-6, -197, -22})),
        baf::ControlPoint(std::vector<double>({-47, -96, -33})),
        baf::ControlPoint(std::vector<double>({25, 32, 95})),
        baf::ControlPoint(std::vector<double>({-6, 102, 8})),
        baf::ControlPoint(std::vector<double>({74, -177, 8})),
        baf::ControlPoint(std::vector<double>({34, -75, 147})),
        baf::ControlPoint(std::vector<double>({86, 97, 105})),
        baf::ControlPoint(std::vector<double>({54, 142, 8})),
        baf::ControlPoint(std::vector<double>({124, -157, 8})),
        baf::ControlPoint(std::vector<double>({198, -31, 63})),
        baf::ControlPoint(std::vector<double>({64, 31, 154})),
        baf::ControlPoint(std::vector<double>({144, 102, 8})),
        baf::ControlPoint(std::vector<double>({204, -197, -22})),
        baf::ControlPoint(std::vector<double>({234, -77, 8})),
        baf::ControlPoint(std::vector<double>({214, -7, 8})),
        baf::ControlPoint(std::vector<double>({239, 102, -22}))
    };
    std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector_ptr =
        {std::make_shared<baf::KnotVector>(knot_vector[0]), std::make_shared<baf::KnotVector>(knot_vector[1])};
    b_spline_ = std::make_unique<spl::BSpline<2>>(knot_vector_ptr, degree, control_points);
  }

 protected:
  std::unique_ptr<spl::Spline<2>> b_spline_;
};

TEST_F(ABSplineSurface, CloseToCenter) { // NOLINT
  std::array<double, 2>
      parametric_coordinates = spl::Projection<2>::ProjectionOnSurface({120, 10, 100}, std::move(b_spline_));
  ASSERT_THAT(parametric_coordinates[0], DoubleNear(0.55852, 0.000001));
  ASSERT_THAT(parametric_coordinates[1], DoubleNear(0.8614466, 0.000001));
}
