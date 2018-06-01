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

class A2DNurbs : public Test {
 public:
  A2DNurbs() {

    std::array<baf::KnotVector, 2> knot_vector = {baf::KnotVector({0, 0, 0, 1, 1, 1}),
                                                  baf::KnotVector({0, 0, 0, 1, 1, 1})};
    std::array<int, 2> degree = {2, 2};
    std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1, 2, 1};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({3.0, 0.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 0.5})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.5, 3.5})),
        baf::ControlPoint(std::vector<double>({5.0, 2.0}))
    };
    nurbs_ = std::make_unique<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::unique_ptr<spl::NURBS<2>> nurbs_;
};

TEST_F(A2DNurbs, Returns1_6For0_4And0_6AndDim0) {
  ASSERT_THAT(nurbs_->Evaluate({0.4, 0.6}, {0})[0], DoubleNear(1.62074, 0.00001));
}

TEST_F(A2DNurbs, Returns1_9For0_4And0_6AndDim1) {
  ASSERT_THAT(nurbs_->Evaluate({0.4, 0.6}, {1})[0], DoubleNear(1.88267, 0.00001));
}

TEST_F(A2DNurbs, Returns2_5For0_5And1_0AndDim0) {
  ASSERT_THAT(nurbs_->Evaluate({0.5, 1.0}, {0})[0], DoubleNear(2.5, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(A2DNurbs, Returns3_0For0_5And1_0AndDim1) {
  ASSERT_THAT(nurbs_->Evaluate({0.5, 1.0}, {1})[0], DoubleNear(3.0, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(A2DNurbs, Returns4_2For0_9And1_0AndDim0) {
  ASSERT_THAT(nurbs_->Evaluate({0.9, 1.0}, {0})[0], DoubleNear(4.19492, 0.00001));
}

TEST_F(A2DNurbs, Returns2_5For0_9And1_0AndDim1) {
  ASSERT_THAT(nurbs_->Evaluate({0.9, 1.0}, {1})[0], DoubleNear(2.45763, 0.00001));
}

TEST_F(A2DNurbs, EvaluatesMultipleValues) {
  ASSERT_THAT(nurbs_->Evaluate({1.0, 0.0}, {0, 1})[0], DoubleEq(3.0));
  ASSERT_THAT(nurbs_->Evaluate({1.0, 0.0}, {0, 1})[1], DoubleEq(0.0));
}

TEST_F(A2DNurbs, Returns10_0For0_0And1_0ForDerivative1And0AndDim0) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({0.0, 1.0}, {0}, {1, 0})[0], DoubleEq(10.0));
}

TEST_F(A2DNurbs, Returns6_0For0_0And1_0ForDerivative1And0AndDim1) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({0.0, 1.0}, {1}, {1, 0})[0], DoubleEq(6.0));
}

TEST_F(A2DNurbs, Returns2_0For0_0And1_0ForDerivative0And1AndDim0) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({0.0, 1.0}, {0}, {0, 1})[0], DoubleEq(2.0));
}

TEST_F(A2DNurbs, Returns3_0For0_0And1_0ForDerivative0And1AndDim1) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({0.0, 1.0}, {1}, {0, 1})[0], DoubleEq(3.0));
}

TEST_F(A2DNurbs, Returns4_2For0_4And0_6ForDerivative1And0AndDim0) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({0.4, 0.6}, {0}, {1, 0})[0], DoubleNear(4.15298, 0.000001));
}

TEST_F(A2DNurbs, Returns0_8For0_4And0_6ForDerivative1And0AndDim1) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({0.4, 0.6}, {1}, {1, 0})[0], DoubleNear(0.792032, 0.000001));
}

TEST_F(A2DNurbs, Returns1_4For0_4And0_6ForDerivative0And1AndDim0) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({0.4, 0.6}, {0}, {0, 1})[0], DoubleNear(1.40046, 0.00001));
}

TEST_F(A2DNurbs, Returns3_1For0_4And0_6ForDerivative0And1AndDim1) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({0.4, 0.6}, {1}, {0, 1})[0], DoubleNear(3.13402, 0.00001));
}

class A2DNurbsWithAllWeights1 : public Test {
 public:
  A2DNurbsWithAllWeights1() {

    std::array<baf::KnotVector, 2> knot_vector = {baf::KnotVector({0, 0, 0, 1, 1, 1}),
                                                  baf::KnotVector({0, 0, 0, 1, 1, 1})};
    std::array<int, 2> degree = {2, 2};
    std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({3.0, 0.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 0.5})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.5, 3.5})),
        baf::ControlPoint(std::vector<double>({5.0, 2.0}))
    };
    nurbs_ = std::make_unique<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
    bspline_ = std::make_unique<spl::BSpline<2>>(knot_vector, degree, control_points);
  }

 protected:
  std::unique_ptr<spl::NURBS<2>> nurbs_;
  std::unique_ptr<spl::BSpline<2>> bspline_;
};

TEST_F(A2DNurbsWithAllWeights1, ReturnsSameDerivativeAs2DBSplineFor0_5And0_5AndDerivatives1And1) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({0.5, 0.5}, {0}, {1, 1})[0],
              DoubleEq(bspline_->EvaluateDerivative({0.5, 0.5}, {0}, {1, 1})[0]));
}

TEST_F(A2DNurbsWithAllWeights1, ReturnsSameDerivativeAs2DBSplineFor0_0And0_7AndDerivatives1And1) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({0.0, 0.7}, {0}, {1, 1})[0],
              DoubleEq(bspline_->EvaluateDerivative({0.0, 0.7}, {0}, {1, 1})[0]));
}

TEST_F(A2DNurbsWithAllWeights1, ReturnsSameDerivativeAs2DBSplineFor0_0And0_7AndDerivatives2And1) {
  ASSERT_THAT(nurbs_->EvaluateDerivative({0.0, 0.7}, {0}, {2, 1})[0],
              DoubleNear(bspline_->EvaluateDerivative({0.0, 0.7}, {0}, {2, 1})[0], 0.000001));
}
