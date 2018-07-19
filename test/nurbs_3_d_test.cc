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
using testing::DoubleEq;
using testing::DoubleNear;

class A3DNurbsWithAllWeights1 : public Test {
 public:
  A3DNurbsWithAllWeights1() {
    std::array<baf::KnotVector, 3> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}),
         baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}}),
         baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})};
    std::array<int, 3> degree = {2, 1, 1};
    std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({3.0, 0.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 0.5})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 0.5})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({5.0, 2.0}))
    };
    std::shared_ptr<std::array<baf::KnotVector, 3>>
        knot_vector_ptr = std::make_shared<std::array<baf::KnotVector, 3>>(knot_vector);
    nurbs_ = std::make_unique<spl::NURBS<3>>(knot_vector_ptr, degree, control_points, weights);
    bspline_ = std::make_unique<spl::BSpline<3>>(knot_vector_ptr, degree, control_points);
  }

 protected:
  std::unique_ptr<spl::NURBS<3>> nurbs_;
  std::unique_ptr<spl::BSpline<3>> bspline_;
};

TEST_F(A3DNurbsWithAllWeights1, ReturnsSameDerivativeAs3DBSplineFor0_5And0_5And0_5AndDerivatives1And1And0) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}, ParamCoord{0.5}}, {0}, {1, 1, 0})[0],
              DoubleEq(bspline_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}, ParamCoord{0.5}},
                                                    {0},
                                                    {1, 1, 0})[0]));
}

TEST_F(A3DNurbsWithAllWeights1, ReturnsSameDerivativeAs3DBSplineFor0_5And0_8And0_1AndDerivatives1And1And1) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.8}, ParamCoord{0.1}}, {0}, {1, 1, 1})[0],
              DoubleEq(bspline_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.8}, ParamCoord{0.1}},
                                                    {0},
                                                    {1, 1, 1})[0]));
}

TEST_F(A3DNurbsWithAllWeights1, ReturnsSameDerivativeAs3DBSplineFor0_5And0_8And0_1AndDerivatives1And2And1) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.8}, ParamCoord{0.1}}, {0}, {1, 2, 1})[0],
              DoubleNear(bspline_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.8}, ParamCoord{0.1}},
                                                      {0},
                                                      {1, 2, 1})[0],
                         util::NumericSettings<double>::kEpsilon()));
}
