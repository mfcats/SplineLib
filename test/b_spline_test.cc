/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <array>
#include <numeric>

#include "gmock/gmock.h"

#include "b_spline.h"
#include "b_spline_generator.h"
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
  ABSpline() : degree_{2} {
    std::array<baf::KnotVector, 1> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2}, ParamCoord{3},
                          ParamCoord{4}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}})};
    control_points_ = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.5, 1.5})),
        baf::ControlPoint(std::vector<double>({2.0, 1.3})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.5})),
        baf::ControlPoint(std::vector<double>({4.0, 0.0}))
    };
    knot_vector_ = std::make_shared<std::array<baf::KnotVector, 1>>(knot_vector);
    b_spline = std::make_unique<spl::BSpline<1>>(knot_vector_, degree_, control_points_);
  }

 protected:
  std::unique_ptr<spl::BSpline<1>> b_spline;
  std::shared_ptr<std::array<baf::KnotVector, 1>> knot_vector_;
  std::array<int, 1> degree_;
  std::vector<baf::ControlPoint> control_points_;
};

TEST_F(ABSpline, Returns0_0For0AndDim0) { // NOLINT
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{0.0}}, {0})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns0_0For0AndDim1) { // NOLINT
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{0.0}}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns4_0For5AndDim0) { // NOLINT
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{5.0}}, {0})[0], DoubleEq(4.0));
}

TEST_F(ABSpline, Returns0_0For5AndDim1) { // NOLINT
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{5.0}}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns1_5For2_5AndDim0) { // NOLINT
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{2.5}}, {0})[0], DoubleEq(1.5));
}

TEST_F(ABSpline, Returns0_0For0_0Dim0AndDer1) { // NOLINT
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{0.0}}, {0}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns1_0For0_0Dim1AndDer1) { // NOLINT
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{0.0}}, {1}, {1})[0], DoubleEq(2.0));
}

TEST_F(ABSpline, Returns12_0For5_0Dim0AndDer1) { // NOLINT
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{5.0}}, {0}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns0_325For2_25Dim1AndDer1) { // NOLINT
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{2.25}}, {1}, {1})[0], DoubleEq(0.325));
}

TEST_F(ABSpline, CanBeConstructedWithAPhysicalAndAParameterSpace) { // NOLINT
  spl::ParameterSpace<1> parameter_space = spl::ParameterSpace<1>({(*knot_vector_)[0]}, {degree_[0]});
  spl::PhysicalSpace<1> physical_space = spl::PhysicalSpace<1>(control_points_, {8});
  b_spline = std::make_unique<spl::BSpline<1>>(parameter_space, physical_space);
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{5.0}}, {0})[0], DoubleEq(4.0));
}

TEST_F(ABSpline, ThrowsExceptionForDifferentControlPointDimensions) { // NOLINT
  control_points_.emplace_back(std::vector<double>({0.0}));
  ASSERT_THROW(std::make_unique<spl::BSpline<1>>(knot_vector_, degree_, control_points_), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForNotFittingKnotVectorAndControlPoints) { // NOLINT
  control_points_.emplace_back(std::vector<double>({0.0, 1.0}));
  ASSERT_THROW(std::make_unique<spl::BSpline<1>>(knot_vector_, degree_, control_points_), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForEvaluationAt6_0) { // NOLINT
  ASSERT_THROW(b_spline->Evaluate({ParamCoord{6.0}}, {0}), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForEvaluationAtMinus1_0) { // NOLINT
  ASSERT_THROW(b_spline->Evaluate({ParamCoord{-1.0}}, {0}), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForDerivativeEvaluationAt6_0) { // NOLINT
  ASSERT_THROW(b_spline->EvaluateDerivative({ParamCoord{6.0}}, {0}, {1}), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForDerivativeEvaluationAtMinus1_0) { // NOLINT
  ASSERT_THROW(b_spline->EvaluateDerivative({ParamCoord{-1.0}}, {0}, {1}), std::runtime_error);
}

TEST_F(ABSpline, ReturnsCorrectNumberOfElements) { // NOLINT
  ASSERT_THAT(b_spline->GetElementList().size(), 5);
}

TEST_F(ABSpline, Returns1DElements) { // NOLINT
  for (auto &element : b_spline->GetElementList()) {
    ASSERT_THAT(element.GetDimension(), 1);
  }
}

TEST_F(ABSpline, ReturnsElementsWith2Nodes) { // NOLINT
  for (auto &element : b_spline->GetElementList()) {
    ASSERT_THAT(element.GetNumberOfNodes(), 2);
  }
}

TEST_F(ABSpline, ReturnsElementsWithCorrectNodes) { // NOLINT
  auto element_list = b_spline->GetElementList();
  for (auto element = 0u; element < element_list.size(); element++) {
    ASSERT_THAT(element_list[element].GetNode(0), element);
    ASSERT_THAT(element_list[element].GetNode(1), element + 1);
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

TEST_F(ABSpline, ReturnsCorrectNonZeroElementBasisFunctionsFor1PointIntegrationRule) { // NOLINT
  auto values = b_spline->EvaluateAllElementNonZeroBasisFunctions(1, itg::IntegrationRule<1>(
      itg::OnePointGaussLegendre<1>()));
  ASSERT_THAT(values.size(), 1);
  ASSERT_THAT(values[0].GetNumberOfNonZeroBasisFunctions(), 3);
  ASSERT_THAT(values[0].GetBasisFunctionValue(0), DoubleEq(0.125));
  ASSERT_THAT(values[0].GetBasisFunctionValue(1), DoubleEq(0.75));
  ASSERT_THAT(values[0].GetBasisFunctionValue(2), DoubleEq(0.125));
}

TEST_F(AnIntegrationRule, LeadsToCorrectNumberOfNonZeroElementBasisFunctions) { // NOLINT
  for (int rule = 1; rule <= 5; rule++) {
    for (int element = 0; element < 5; element++) {
      auto values = b_spline->EvaluateAllElementNonZeroBasisFunctions(element, rules_[rule - 1]);
      ASSERT_THAT(values.size(), rule);
      for (int point = 0; point < rule; point++) {
        ASSERT_THAT(values[point].GetNumberOfNonZeroBasisFunctions(), 3);
        std::vector<double> non_zero_basis_functions = values[point].GetNonZeroBasisFunctions();
        ASSERT_THAT(std::accumulate(non_zero_basis_functions.cbegin(),
                                    non_zero_basis_functions.cend(),
                                    0.0,
                                    std::plus<>()), DoubleEq(1.0));
      }
    }
  }
}

TEST_F(ABSpline, ReturnsCorrectNonZeroElementBasisFunctionDerivativesFor1PointIntegrationRule) { // NOLINT
  auto values =
      b_spline->EvaluateAllElementNonZeroBasisFunctionDerivatives(1, itg::IntegrationRule<1>(
          itg::OnePointGaussLegendre<1>()));
  ASSERT_THAT(values.size(), 1);
  ASSERT_THAT(values[0].GetNumberOfNonZeroBasisFunctions(), 3);
  ASSERT_THAT(values[0].GetBasisFunctionValue(0), DoubleEq(-2.0 / 3.0));
  ASSERT_THAT(values[0].GetBasisFunctionValue(1), DoubleEq(0));
  ASSERT_THAT(values[0].GetBasisFunctionValue(2), DoubleEq(2.0 / 3.0));
}

TEST_F(AnIntegrationRule, LeadsToCorrectNumberOfNonZeroElementBasisFunctionDerivatives) { // NOLINT
  for (int rule = 1; rule <= 5; rule++) {
    for (int element = 0; element < 5; element++) {
      auto values = b_spline->EvaluateAllElementNonZeroBasisFunctionDerivatives(element, rules_[rule - 1]);
      ASSERT_THAT(values.size(), rule);
      for (int point = 0; point < rule; point++) {
        ASSERT_THAT(values[point].GetNumberOfNonZeroBasisFunctions(), 3);
        std::vector<double> non_zero_basis_functions = values[point].GetNonZeroBasisFunctions();
        ASSERT_THAT(std::accumulate(non_zero_basis_functions.cbegin(),
                                    non_zero_basis_functions.cend(),
                                    0.0,
                                    std::plus<>()),
                    DoubleNear(0.0, util::NumericSettings<double>::kEpsilon()));
      }
    }
  }
}

TEST_F(ABSpline, ReturnsCorrectJacobianDeterminant) { // NOLINT
  ASSERT_THAT(b_spline->JacobianDeterminant(1, 1, itg::ThreePointGaussLegendre<1>()), DoubleEq(0.375));
}

class ABSplineWithSplineGenerator : public Test {
 public:
  ABSplineWithSplineGenerator() {
    std::array<baf::KnotVector, 1> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2}, ParamCoord{3},
                          ParamCoord{4}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}})};
    degree_ = {2};
    control_points_ = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.5, 1.5})),
        baf::ControlPoint(std::vector<double>({2.0, 1.3})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.5})),
        baf::ControlPoint(std::vector<double>({4.0, 0.0}))
    };
    knot_vector_[0] = {std::make_shared<baf::KnotVector>(knot_vector[0])};
    spl::BSplineGenerator<1> b_spline_generator(knot_vector_, degree_, control_points_);
    b_spline = std::make_unique<spl::BSpline<1>>(b_spline_generator);
  }

 protected:
  std::unique_ptr<spl::BSpline<1>> b_spline;
  std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_;
  std::array<int, 1> degree_;
  std::vector<baf::ControlPoint> control_points_;
};

TEST_F(ABSplineWithSplineGenerator, Returns0_0For0AndDim0) { // NOLINT
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{0.0}}, {0})[0], DoubleEq(0.0));
}
