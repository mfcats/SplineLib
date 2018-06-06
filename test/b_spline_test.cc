/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <numeric>

#include "gmock/gmock.h"

#include "b_spline.h"
#include "knot_vector.h"
#include "one_point_gauss_legendre.h"
#include "two_point_gauss_legendre.h"
#include "three_point_gauss_legendre.h"
#include "four_point_gauss_legendre.h"
#include "five_point_gauss_legendre.h"
#include "numeric_settings.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class ABSpline : public Test {
 public:
  ABSpline() {
    std::array<baf::KnotVector, 1> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2}, ParamCoord{3},
                          ParamCoord{4}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}})};
    std::array<int, 1> degree = {2};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.5, 1.5})),
        baf::ControlPoint(std::vector<double>({2.0, 1.3})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.5})),
        baf::ControlPoint(std::vector<double>({4.0, 0.0}))
    };
    b_spline = std::make_unique<spl::BSpline<1>>(knot_vector, degree, control_points);
  }

 protected:
  std::unique_ptr<spl::BSpline<1>> b_spline;
};

TEST_F(ABSpline, Returns0_0For0AndDim0) {
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{0.0}}, {0})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns0_0For0AndDim1) {
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{0.0}}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns4_0For5AndDim0) {
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{5.0}}, {0})[0], DoubleEq(4.0));
}

TEST_F(ABSpline, Returns0_0For5AndDim1) {
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{5.0}}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns1_5For2_5AndDim0) {
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{2.5}}, {0})[0], DoubleEq(1.5));
}

TEST_F(ABSpline, Returns0_0For0_0Dim0AndDer1) {
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{0.0}}, {0}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns1_0For0_0Dim1AndDer1) {
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{0.0}}, {1}, {1})[0], DoubleEq(2.0));
}

TEST_F(ABSpline, Returns12_0For5_0Dim0AndDer1) {
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{5.0}}, {0}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns0_325For2_25Dim1AndDer1) {
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{2.25}}, {1}, {1})[0], DoubleEq(0.325));
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

class AnIntegrationRule : public ABSpline {
 public:
  AnIntegrationRule() {
    rules_.emplace_back(itg::OnePointGaussLegendre<1>());
    rules_.emplace_back(itg::TwoPointGaussLegendre<1>());
    rules_.emplace_back(itg::ThreePointGaussLegendre<1>());
    rules_.emplace_back(itg::FourPointGaussLegendre<1>());
    rules_.emplace_back(itg::FivePointGaussLegendre<1>());
  }

 protected:
  std::vector<itg::IntegrationRule<1>> rules_;
};

TEST_F(ABSpline, ReturnsCorrectNonZeroElementBasisFunctionsFor1PointIntegrationRule) {
  auto values = b_spline->EvaluateAllElementNonZeroBasisFunctions(1, itg::IntegrationRule<1>(
      itg::OnePointGaussLegendre<1>()));
  ASSERT_THAT(values.size(), 1);
  ASSERT_THAT(values[0].NumberOfNonZeroBasisFunctions(), 3);
  ASSERT_THAT(values[0].GetBasisFunctionValue(0), DoubleEq(0.125));
  ASSERT_THAT(values[0].GetBasisFunctionValue(1), DoubleEq(0.75));
  ASSERT_THAT(values[0].GetBasisFunctionValue(2), DoubleEq(0.125));
}

TEST_F(AnIntegrationRule, LeadsToCorrectNumberOfNonZeroElementBasisFunctions) {
  for (int rule = 1; rule <= 5; rule++) {
    for (int element = 0; element < 5; element++) {
      auto values = b_spline->EvaluateAllElementNonZeroBasisFunctions(element, rules_[rule - 1]);
      ASSERT_THAT(values.size(), rule);
      for (int point = 0; point < rule; point++) {
        ASSERT_THAT(values[point].NumberOfNonZeroBasisFunctions(), 3);
        std::vector<double> non_zero_basis_functions = values[point].non_zero_basis_functions();
        ASSERT_THAT(std::accumulate(non_zero_basis_functions.cbegin(),
                                    non_zero_basis_functions.cend(),
                                    0.0,
                                    std::plus<double>()), DoubleEq(1.0));
      }
    }
  }
}

TEST_F(ABSpline, ReturnsCorrectNonZeroElementBasisFunctionDerivativesFor1PointIntegrationRule) {
  auto values =
      b_spline->EvaluateAllElementNonZeroBasisFunctionDerivatives(1, itg::IntegrationRule<1>(
          itg::OnePointGaussLegendre<1>()));
  ASSERT_THAT(values.size(), 1);
  ASSERT_THAT(values[0].NumberOfNonZeroBasisFunctions(), 3);
  ASSERT_THAT(values[0].GetBasisFunctionValue(0), DoubleEq(-2.0 / 3.0));
  ASSERT_THAT(values[0].GetBasisFunctionValue(1), DoubleEq(0));
  ASSERT_THAT(values[0].GetBasisFunctionValue(2), DoubleEq(2.0 / 3.0));
}

TEST_F(AnIntegrationRule, LeadsToCorrectNumberOfNonZeroElementBasisFunctionDerivatives) {
  for (int rule = 1; rule <= 5; rule++) {
    for (int element = 0; element < 5; element++) {
      auto values = b_spline->EvaluateAllElementNonZeroBasisFunctionDerivatives(element, rules_[rule - 1]);
      ASSERT_THAT(values.size(), rule);
      for (int point = 0; point < rule; point++) {
        ASSERT_THAT(values[point].NumberOfNonZeroBasisFunctions(), 3);
        std::vector<double> non_zero_basis_functions = values[point].non_zero_basis_functions();
        ASSERT_THAT(std::accumulate(non_zero_basis_functions.cbegin(),
                                    non_zero_basis_functions.cend(),
                                    0.0,
                                    std::plus<double>()),
                    DoubleNear(0.0, util::NumericSettings<double>::kEpsilon()));
      }
    }
  }
}

TEST_F(ABSpline, ReturnsCorrectJacobianDeterminant) {
  ASSERT_THAT(b_spline->JacobianDeterminant(1, 1, itg::ThreePointGaussLegendre<1>()), DoubleEq(0.375));
}
