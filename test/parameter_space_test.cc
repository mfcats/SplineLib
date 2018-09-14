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

#include "parameter_space.h"
#include "five_point_gauss_legendre.h"
#include "four_point_gauss_legendre.h"
#include "numeric_settings.h"
#include "one_point_gauss_legendre.h"
#include "two_point_gauss_legendre.h"
#include "three_point_gauss_legendre.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class A1DParameterSpace : public Test {
 public:
  A1DParameterSpace() :
      degree_{Degree{2}},
      knot_vector_{
          std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                             ParamCoord{2}, ParamCoord{3}, ParamCoord{4}, ParamCoord{4},
                                                             ParamCoord{5}, ParamCoord{5}, ParamCoord{5}}))},
      parameter_space(spl::ParameterSpace<1>(knot_vector_, degree_)) {}

 protected:
  std::array<Degree, 1> degree_;
  std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_;
  spl::ParameterSpace<1> parameter_space;
};

TEST_F(A1DParameterSpace, returnsCorrectDegree) {  // NOLINT
  ASSERT_THAT(parameter_space.GetDegree(0), Degree{2});
}

TEST_F(A1DParameterSpace, returns3_0ForFifthKnot) {  // NOLINT
  ASSERT_THAT((*parameter_space.GetKnotVector(0))[5].get(), DoubleEq(3.0));
}

TEST_F(A1DParameterSpace, returnsCorrectBasisFunctionValuesForParamCoord0_5) {  // NOLINT
  std::vector<double> values = {0.25, 0.625, 0.125, 0, 0, 0, 0, 0};
  for (int i = 0; i < 8; ++i) {
    ASSERT_THAT(parameter_space.GetBasisFunctions({i}, {ParamCoord(0.5)}), DoubleEq(values[i]));
  }
}

TEST_F(A1DParameterSpace, returnsCorrectBasisFunctionDerivativeValuesForParamCoord0_5AndDerivative1) {  // NOLINT
  std::vector<double> values = {-1, 0.5, 0.5, 0, 0, 0, 0, 0};
  for (int i = 0; i < 8; ++i) {
    ASSERT_THAT(parameter_space.GetBasisFunctionDerivatives({i}, {ParamCoord(0.5)}, {1}), DoubleEq(values[i]));
  }
}

TEST_F(A1DParameterSpace, ReturnsCorrectNumberOfElements) { // NOLINT
  ASSERT_THAT(parameter_space.GetElementList(0).size(), 5);
}

TEST_F(A1DParameterSpace, Returns1DElements) { // NOLINT
  for (auto &element : parameter_space.GetElementList(0)) {
    ASSERT_THAT(element.GetDimension(), 1);
  }
}

TEST_F(A1DParameterSpace, ReturnsElementsWith2Nodes) { // NOLINT
  for (auto &element : parameter_space.GetElementList(0)) {
    ASSERT_THAT(element.GetNumberOfNodes(), 2);
  }
}

TEST_F(A1DParameterSpace, ReturnsElementsWithCorrectNodes) { // NOLINT
  auto element_list = parameter_space.GetElementList(0);
  for (auto element = 0u; element < element_list.size(); element++) {
    ASSERT_THAT(element_list[element].GetNode(0).get(), element);
    ASSERT_THAT(element_list[element].GetNode(1).get(), element + 1);
  }
}

TEST_F(A1DParameterSpace, ReturnsCorrectNonZeroElementBasisFunctionsFor1PointIntegrationRule) { // NOLINT
  auto values = parameter_space.EvaluateAllElementNonZeroBasisFunctions(0, 1, itg::IntegrationRule<1>(
      itg::OnePointGaussLegendre<1>()));
  ASSERT_THAT(values.size(), 1);
  ASSERT_THAT(values[0].GetNumberOfNonZeroBasisFunctions(), 3);
  ASSERT_THAT(values[0].GetBasisFunctionValue(0), DoubleEq(0.125));
  ASSERT_THAT(values[0].GetBasisFunctionValue(1), DoubleEq(0.75));
  ASSERT_THAT(values[0].GetBasisFunctionValue(2), DoubleEq(0.125));
}

TEST_F(A1DParameterSpace, ReturnsCorrectNonZeroElementBasisFunctionDerivativesFor1PointIntegrationRule) { // NOLINT
  auto values = parameter_space.EvaluateAllElementNonZeroBasisFunctionDerivatives(0, 1, itg::IntegrationRule<1>(
      itg::OnePointGaussLegendre<1>()));
  ASSERT_THAT(values.size(), 1);
  ASSERT_THAT(values[0].GetNumberOfNonZeroBasisFunctions(), 3);
  ASSERT_THAT(values[0].GetBasisFunctionValue(0), DoubleEq(-0.5));
  ASSERT_THAT(values[0].GetBasisFunctionValue(1), DoubleEq(0));
  ASSERT_THAT(values[0].GetBasisFunctionValue(2), DoubleEq(0.5));
}

TEST_F(A1DParameterSpace, TransformsReferenceElementCoordinateToParametricCoordinate) { // NOLINT
  ASSERT_THAT(parameter_space.ReferenceSpace2ParameterSpace(ParamCoord{3}, ParamCoord{2}, 0.5).get(), DoubleEq(2.75));
}

TEST_F(A1DParameterSpace, Returns1AsFirstNonZeroBasisFunctionsForParamCoord1_5) { // NOLINT
  ASSERT_THAT(parameter_space.GetArrayOfFirstNonZeroBasisFunctions({ParamCoord(1.5)})[0], 1);
}

TEST_F(A1DParameterSpace, Returns5AsFirstNonZeroBasisFunctionsForParamCoord4) { // NOLINT
  ASSERT_THAT(parameter_space.GetArrayOfFirstNonZeroBasisFunctions({ParamCoord(4.0)})[0], 5);
}

TEST_F(A1DParameterSpace, ThrowsIfParametricCoordinateOutsideKnotVectorRange) { // NOLINT
  ASSERT_THROW(parameter_space.ThrowIfParametricCoordinateOutsideKnotVectorRange({ParamCoord(5.1)}), std::range_error);
}

class AnIntegrationRule : public A1DParameterSpace {
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

TEST_F(AnIntegrationRule, LeadsToCorrectNumberOfNonZeroElementBasisFunctions) { // NOLINT
  for (int rule = 1; rule <= 5; rule++) {
    for (int element = 0; element < 5; element++) {
      auto values = parameter_space.EvaluateAllElementNonZeroBasisFunctions(0, element, rules_[rule - 1]);
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

TEST_F(AnIntegrationRule, LeadsToCorrectNumberOfNonZeroElementBasisFunctionDerivatives) { // NOLINT
  for (int rule = 1; rule <= 5; rule++) {
    for (int element = 0; element < 5; element++) {
      auto values = parameter_space.EvaluateAllElementNonZeroBasisFunctionDerivatives(0, element, rules_[rule - 1]);
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

class A2DParameterSpace : public Test {
 public:
  A2DParameterSpace() : degree_{Degree{2}, Degree{1}},
                        knot_vector_{std::make_shared<baf::KnotVector>(baf::KnotVector{ParamCoord{0}, ParamCoord{0},
                                                                                       ParamCoord{0}, ParamCoord{1},
                                                                                       ParamCoord{1}, ParamCoord{1}}),
                                     std::make_shared<baf::KnotVector>(baf::KnotVector{ParamCoord{0}, ParamCoord{0},
                                                                                       ParamCoord{1}, ParamCoord{2},
                                                                                       ParamCoord{3}, ParamCoord{3}})},
                        parameter_space(knot_vector_, degree_) {}

 protected:
  std::array<Degree, 2> degree_;
  std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector_;
  spl::ParameterSpace<2> parameter_space;
};

TEST_F(A2DParameterSpace, returnsDegree2ForFirstDimension) {  // NOLINT
  ASSERT_THAT(parameter_space.GetDegree(0), Degree{2});
}

TEST_F(A2DParameterSpace, returnsDegree1ForSecondDimension) {  // NOLINT
  ASSERT_THAT(parameter_space.GetDegree(1), Degree{1});
}

TEST_F(A2DParameterSpace, returns1_0ForFourthKnotOfFirstKnotVector) {  // NOLINT
  ASSERT_THAT((*parameter_space.GetKnotVector(0))[3].get(), DoubleEq(1.0));
}

TEST_F(A2DParameterSpace, returns2_0ForFourthKnotOfSecondKnotVector) {  // NOLINT
  ASSERT_THAT((*parameter_space.GetKnotVector(1))[3].get(), DoubleEq(2.0));
}

TEST_F(A2DParameterSpace, returnsCorrectBasisFunctionValue) {  // NOLINT
  ASSERT_THAT(parameter_space.GetBasisFunctions({1, 1}, {ParamCoord(0.5), ParamCoord(0.25)}), DoubleEq(0.125));
}

TEST_F(A2DParameterSpace, returnsCorrectBasisFunctionDerivativeValue) {  // NOLINT
  ASSERT_THAT(parameter_space.GetBasisFunctionDerivatives({1, 1}, {ParamCoord(0.5), ParamCoord(0.25)}, {0, 1}),
              DoubleEq(0.5));
}

class A3DParameterSpace : public Test {
 public:
  A3DParameterSpace() : degree_{Degree{2}, Degree{0}, Degree{1}},
                        knot_vector_{
                            std::make_shared<baf::KnotVector>(baf::KnotVector(
                                {std::vector<ParamCoord>({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                          ParamCoord{1}, ParamCoord{1}})})),
                            std::make_shared<baf::KnotVector>(baf::KnotVector(
                                {std::vector<ParamCoord>({ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6},
                                                          ParamCoord{0.9}})})),
                            std::make_shared<baf::KnotVector>(baf::KnotVector(
                                {std::vector<ParamCoord>({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2},
                                                          ParamCoord{3}, ParamCoord{3}})}))},
                        parameter_space(knot_vector_, degree_) {}

 protected:
  std::array<Degree, 3> degree_;
  std::array<std::shared_ptr<baf::KnotVector>, 3> knot_vector_;
  spl::ParameterSpace<3> parameter_space;
};

TEST_F(A3DParameterSpace, returnsDegree2ForFirstDimension) {  // NOLINT
  ASSERT_THAT(parameter_space.GetDegree(0), Degree{2});
}

TEST_F(A3DParameterSpace, returnsDegree0ForSecondDimension) {  // NOLINT
  ASSERT_THAT(parameter_space.GetDegree(1), Degree{0});
}

TEST_F(A3DParameterSpace, returnsDegree1ForThirdDimension) {  // NOLINT
  ASSERT_THAT(parameter_space.GetDegree(2), Degree{1});
}

TEST_F(A3DParameterSpace, returns1_0ForFourthKnotOfFirstKnotVector) {  // NOLINT
  ASSERT_THAT((*parameter_space.GetKnotVector(0))[3].get(), DoubleEq(1.0));
}

TEST_F(A3DParameterSpace, returns0_3ForSecondKnotOfSecondKnotVector) {  // NOLINT
  ASSERT_THAT((*parameter_space.GetKnotVector(1))[1].get(), DoubleEq(0.3));
}

TEST_F(A3DParameterSpace, returns2_0ForFourthKnotOfThirdKnotVector) {  // NOLINT
  ASSERT_THAT((*parameter_space.GetKnotVector(2))[3].get(), DoubleEq(2.0));
}

TEST_F(A3DParameterSpace, returnsCorrectBasisFunctionValue) {  // NOLINT
  ASSERT_THAT(parameter_space.GetBasisFunctions({1, 1, 1}, {ParamCoord(0.5), ParamCoord(0.5), ParamCoord(0.25)}),
              DoubleEq(0.125));
}

TEST_F(A3DParameterSpace, returnsCorrectBasisFunctionDerivativeValue) {  // NOLINT
  ASSERT_THAT(parameter_space.GetBasisFunctionDerivatives({1, 1, 1},
                                                          {ParamCoord(0.5), ParamCoord(0.5), ParamCoord(0.25)},
                                                          {0, 0, 1}), DoubleEq(0.5));
}

