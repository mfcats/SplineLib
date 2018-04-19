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

class ABSpline : public Test {
 public:
  ABSpline() {
    KnotVector knot_vector = {0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
    std::vector<ControlPoint> control_points = {
        ControlPoint(std::vector<double>({0.0, 0.0})),
        ControlPoint(std::vector<double>({0.0, 1.0})),
        ControlPoint(std::vector<double>({1.0, 1.0})),
        ControlPoint(std::vector<double>({1.5, 1.5})),
        ControlPoint(std::vector<double>({2.0, 1.3})),
        ControlPoint(std::vector<double>({3.0, 2.0})),
        ControlPoint(std::vector<double>({4.0, 1.5})),
        ControlPoint(std::vector<double>({4.0, 0.0}))
    };
    b_spline = std::make_unique<BSpline>(knot_vector, 2, control_points);
  }

 protected:
  std::unique_ptr<BSpline> b_spline;
};

TEST_F(ABSpline, Returns0_0For0AndDim0) {
  ASSERT_THAT(b_spline->Evaluate(0.0, {0})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns0_0For0AndDim1) {
  ASSERT_THAT(b_spline->Evaluate(0.0, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns4_0For5AndDim0) {
  ASSERT_THAT(b_spline->Evaluate(5.0, {0})[0], DoubleEq(4.0));
}

TEST_F(ABSpline, Returns0_0For5AndDim1) {
  ASSERT_THAT(b_spline->Evaluate(5.0, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns1_5For2_5AndDim0) {
  ASSERT_THAT(b_spline->Evaluate(2.5, {0})[0], DoubleEq(1.5));
}

TEST_F(ABSpline, Returns0_0For0_0Dim0AndDer1) {
  ASSERT_THAT(b_spline->EvaluateDerivative(0.0, {0}, 1)[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns1_0For0_0Dim1AndDer1) {
  ASSERT_THAT(b_spline->EvaluateDerivative(0.0, {1}, 1)[0], DoubleEq(2.0));
}

TEST_F(ABSpline, Returns12_0For5_0Dim0AndDer1) {
  ASSERT_THAT(b_spline->EvaluateDerivative(5.0, {0}, 1)[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns0_325For2_25Dim1AndDer1) {
  ASSERT_THAT(b_spline->EvaluateDerivative(2.25, {1}, 1)[0], DoubleEq(0.325));
}

TEST_F(ABSpline, ReturnsCorrectNumberOfElements) {
  ASSERT_THAT(b_spline->GetElementList().size(), 5);
}

TEST_F(ABSpline, Returns1DElements) {
  for (auto &element : b_spline->GetElementList()) {
    ASSERT_THAT(element.dimension(), 1);
  }
}

TEST_F(ABSpline, ReturnsElementsWith2Nodes) {
  for (auto &element : b_spline->GetElementList()) {
    ASSERT_THAT(element.numberOfNodes(), 2);
  }
}

TEST_F(ABSpline, ReturnsElementsWithCorrectNodes) {
  auto element_list = b_spline->GetElementList();
  for (int element = 0; element < element_list.size(); element++) {
    ASSERT_THAT(element_list[element].node(0), element);
    ASSERT_THAT(element_list[element].node(1), element + 1);
  }
}

TEST_F(ABSpline, ReturnsCorrectNumberOfNonZeroElementBasisFunctions) {
  auto values = b_spline->EvaluateAllElementNonZeroBasisFunctions(1, IntegrationRule<1>(1));
  ASSERT_THAT(values.size(), 1);
  ASSERT_THAT(values[0].size(), 3);
  ASSERT_THAT(values[0][0], DoubleEq(0.125));
  ASSERT_THAT(values[0][1], DoubleEq(0.75));
  ASSERT_THAT(values[0][2], DoubleEq(0.125));
}

TEST_F(ABSpline, ReturnsCorrectNumberOfNonZeroElementBasisFunctionsFor2PointIntegrationRule) {
  auto values = b_spline->EvaluateAllElementNonZeroBasisFunctions(1, IntegrationRule<1>(2));
  ASSERT_THAT(values.size(), 2);
  ASSERT_THAT(values[0].size(), 3);
  ASSERT_THAT(values[1].size(), 3);
}
