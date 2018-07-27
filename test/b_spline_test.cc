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
#include <sstream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "b_spline.h"
#include "one_point_gauss_legendre.h"
#include "two_point_gauss_legendre.h"
#include "three_point_gauss_legendre.h"
#include "four_point_gauss_legendre.h"
#include "five_point_gauss_legendre.h"
#include "numeric_settings.h"
#include "b_spline_generator.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;
using ::testing::Return;
using ::testing::Throw;
using ::testing::NiceMock;

class MockParameterSpace : public spl::ParameterSpace<1> {
 public:
  MOCK_CONST_METHOD1(GetDegree, int(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 1>, std::array<ParamCoord, 1>));
  MOCK_CONST_METHOD3(GetBasisFunctionDerivatives, double(std::array<int, 1>, std::array<ParamCoord, 1>, std::array<int, 1>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 1>(std::array<ParamCoord, 1>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParamCoord, 1>));
};

void default_on_call(std::shared_ptr<NiceMock<MockParameterSpace>> parameter_space){
  //Span - degree
  ON_CALL(*parameter_space, ThrowIfParametricCoordinateOutsideKnotVectorRange(std::array<ParamCoord, 1>{ParamCoord{-1.0}}))
      .WillByDefault(Throw(std::range_error("Out of knotvector range")));
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(std::array<ParamCoord, 1>{ParamCoord{0.0}}))
      .WillByDefault(Return(std::array<int, 1>{0}));
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(std::array<ParamCoord, 1>{ParamCoord{2.25}}))
      .WillByDefault(Return(std::array<int, 1>{2}));
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(std::array<ParamCoord, 1>{ParamCoord{2.5}}))
      .WillByDefault(Return(std::array<int, 1>{2}));
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(std::array<ParamCoord, 1>{ParamCoord{5.0}}))
      .WillByDefault(Return(std::array<int, 1>{5}));
  ON_CALL(*parameter_space, ThrowIfParametricCoordinateOutsideKnotVectorRange(std::array<ParamCoord, 1>{ParamCoord{6.0}}))
      .WillByDefault(Throw(std::range_error("Out of knotvector range")));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1> {0}, std::array<ParamCoord, 1>{ParamCoord{0.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1> {2}, std::array<ParamCoord, 1>{ParamCoord{2.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1> {3}, std::array<ParamCoord, 1>{ParamCoord{2.5}}))
      .WillByDefault(Return(0.75));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1> {4}, std::array<ParamCoord, 1>{ParamCoord{2.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1> {5}, std::array<ParamCoord, 1>{ParamCoord{5.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1> {6}, std::array<ParamCoord, 1>{ParamCoord{5.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1> {7}, std::array<ParamCoord, 1>{ParamCoord{5.0}}))
      .WillByDefault(Return(1.0));
  ON_CALL(*parameter_space, GetDegree(0))
      .WillByDefault(Return(2));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1> {0}, std::array<ParamCoord, 1>{ParamCoord{0.0}}, std::array<int, 1>{1}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1> {1}, std::array<ParamCoord, 1>{ParamCoord{0.0}}, std::array<int, 1>{1}))
      .WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1> {2}, std::array<ParamCoord, 1>{ParamCoord{0.0}}, std::array<int, 1>{1}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1> {5}, std::array<ParamCoord, 1>{ParamCoord{5.0}}, std::array<int, 1>{1}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1> {6}, std::array<ParamCoord, 1>{ParamCoord{5.0}}, std::array<int, 1>{1}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1> {7}, std::array<ParamCoord, 1>{ParamCoord{5.0}}, std::array<int, 1>{1}))
      .WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1> {2}, std::array<ParamCoord, 1>{ParamCoord{2.25}}, std::array<int, 1>{1}))
      .WillByDefault(Return(-0.75));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1> {3}, std::array<ParamCoord, 1>{ParamCoord{2.25}}, std::array<int, 1>{1}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1> {4}, std::array<ParamCoord, 1>{ParamCoord{2.25}}, std::array<int, 1>{1}))
      .WillByDefault(Return(0.25));


}
// U = {0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}
class ABSpline : public Test {
 public:
  ABSpline() :
      parameter_space1(std::make_shared<NiceMock<MockParameterSpace>>()) {
    std::array<baf::KnotVector, 1> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2}, ParamCoord{3},
                          ParamCoord{4}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}})};
    std::array<int, 1> degree_ = {2};
    std::vector<baf::ControlPoint> control_points_ = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.5, 1.5})),
        baf::ControlPoint(std::vector<double>({2.0, 1.3})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.5})),
        baf::ControlPoint(std::vector<double>({4.0, 0.0}))
    };
    std::array<std::shared_ptr<baf::KnotVector>, 1>knot_vector_ = {std::make_shared<baf::KnotVector>(knot_vector[0])};

    parameter_space = std::make_shared<spl::ParameterSpace<1>>(spl::ParameterSpace<1>(knot_vector_, degree_));

    std::array<int, 1> number_of_points;
    number_of_points[0] = knot_vector[0].GetNumberOfKnots() - degree_[0] - 1;

    physical_space = std::make_shared<spl::PhysicalSpace<1>>(spl::PhysicalSpace<1>(control_points_, number_of_points));

    spl::BSplineGenerator<1> b_spline_generator1(physical_space, parameter_space);
    spl::BSplineGenerator<1> b_spline_generator(physical_space, parameter_space1);
    b_spline = std::make_unique<spl::BSpline<1>>(b_spline_generator);
    b_spline1 = std::make_unique<spl::BSpline<1>>(b_spline_generator1);
  }

 protected:
  std::unique_ptr<spl::BSpline<1>> b_spline;
  std::unique_ptr<spl::BSpline<1>> b_spline1;
  std::shared_ptr<NiceMock<MockParameterSpace>> parameter_space1;
  std::shared_ptr<spl::ParameterSpace<1>> parameter_space;
  std::shared_ptr<spl::PhysicalSpace<1>> physical_space;
};

TEST_F(ABSpline, Returns0_0For0AndDim0) {
  default_on_call(parameter_space1);
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{0.0}}, {0})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns0_0For0AndDim1) {
  default_on_call(parameter_space1);
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{0.0}}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns4_0For5AndDim0) {
  default_on_call(parameter_space1);
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{5.0}}, {0})[0], DoubleEq(4.0));
}

TEST_F(ABSpline, Returns0_0For5AndDim1) {
  default_on_call(parameter_space1);
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{5.0}}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns1_5For2_5AndDim0) {
  default_on_call(parameter_space1);
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{2.5}}, {0})[0], DoubleEq(1.5));
}

TEST_F(ABSpline, Returns0_0For0_0Dim0AndDer1) {
  default_on_call(parameter_space1);
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{0.0}}, {0}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns1_0For0_0Dim1AndDer1) {
  default_on_call(parameter_space1);
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{0.0}}, {1}, {1})[0], DoubleEq(2.0));
}

TEST_F(ABSpline, Returns12_0For5_0Dim0AndDer1) {
 default_on_call(parameter_space1);
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{5.0}}, {0}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns0_325For2_25Dim1AndDer1) {
  default_on_call(parameter_space1);
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{2.25}}, {1}, {1})[0], DoubleEq(0.325));
}

TEST_F(ABSpline, CanBeConstructedWithAPhysicalAndAParameterSpace) {
  default_on_call(parameter_space1);
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{5.0}}, {0})[0], DoubleEq(4.0));
}

TEST_F(ABSpline, ThrowsExceptionForEvaluationAt6_0) {
  default_on_call(parameter_space1);
  ASSERT_THROW(b_spline->Evaluate({ParamCoord{6.0}}, {0}), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForEvaluationAtMinus1_0) {
  default_on_call(parameter_space1);
  ASSERT_THROW(b_spline->Evaluate({ParamCoord{-1.0}}, {0}), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForDerivativeEvaluationAt6_0) {
  default_on_call(parameter_space1);
  ASSERT_THROW(b_spline->EvaluateDerivative({ParamCoord{6.0}}, {0}, {1}), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForDerivativeEvaluationAtMinus1_0) {
  default_on_call(parameter_space1);
  ASSERT_THROW(b_spline->EvaluateDerivative({ParamCoord{-1.0}}, {0}, {1}), std::runtime_error);
}

/*
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
*/
