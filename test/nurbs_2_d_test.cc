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
using testing::Return;
using testing::DoubleNear;
using ::testing::NiceMock;
using ::testing::Throw;
using ::testing::_;

class MockParameterSpace1 : public spl::ParameterSpace<2> {
 public:
  MOCK_CONST_METHOD1(GetDegree, Degree(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 2>, std::array<ParamCoord, 2>));
  MOCK_CONST_METHOD3(GetBasisFunctionDerivatives,
                     double(std::array<int, 2>, std::array<ParamCoord, 2>, std::array<int, 2>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 2>(std::array<ParamCoord, 2>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParamCoord, 2>));
};

class MockWeightedPhysicalSpace1 : public spl::WeightedPhysicalSpace<2> {
 public:
  MOCK_CONST_METHOD1(GetWeight, double(std::array<int, 2>));
  MOCK_CONST_METHOD1(GetHomogenousControlPoint, baf::ControlPoint(std::array<int, 2>));
  MOCK_CONST_METHOD1(GetControlPoint, baf::ControlPoint(std::array<int, 2>));
};

void set_get_basis_function_nurbs(const std::shared_ptr<NiceMock<MockParameterSpace1>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctions(_, std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{1.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2}, std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{1.0}}))
      .WillByDefault(Return(0.25));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 2}, std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{1.0}}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 2}, std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{1.0}}))
      .WillByDefault(Return(0.25));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0}, std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}}))
      .WillByDefault(Return(0.0576));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 0}, std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}}))
      .WillByDefault(Return(0.0768));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 0}, std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}}))
      .WillByDefault(Return(0.0256));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 1}, std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}}))
      .WillByDefault(Return(0.1728));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 1}, std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}}))
      .WillByDefault(Return(0.2304));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 1}, std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}}))
      .WillByDefault(Return(0.0768));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2}, std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}}))
      .WillByDefault(Return(0.1296));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 2}, std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}}))
      .WillByDefault(Return(0.1728));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 2}, std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}}))
      .WillByDefault(Return(0.0576));
  ON_CALL(*parameter_space, GetBasisFunctions(_, std::array<ParamCoord, 2>{ParamCoord{0.9}, ParamCoord{1.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2}, std::array<ParamCoord, 2>{ParamCoord{0.9}, ParamCoord{1.0}}))
      .WillByDefault(Return(0.01));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 2}, std::array<ParamCoord, 2>{ParamCoord{0.9}, ParamCoord{1.0}}))
      .WillByDefault(Return(0.18));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 2}, std::array<ParamCoord, 2>{ParamCoord{0.9}, ParamCoord{1.0}}))
      .WillByDefault(Return(0.81));
  ON_CALL(*parameter_space, GetBasisFunctions(_, std::array<ParamCoord, 2>{ParamCoord{1.0}, ParamCoord{0.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 0}, std::array<ParamCoord, 2>{ParamCoord{1.0}, ParamCoord{0.0}}))
      .WillByDefault(Return(1));
  ON_CALL(*parameter_space, GetBasisFunctions(_, std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{1.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2}, std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{1.0}}))
      .WillByDefault(Return(1));
}

void set_basis_function_derivative1(const std::shared_ptr<NiceMock<MockParameterSpace1>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(_, _, _)).WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{1.0}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{1.0}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{1.0}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{1.0}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.192));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.064));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.128));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.576));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.192));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.384));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.432));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.144));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.288));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.288));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.384));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.128));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.144));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.192));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.064));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.432));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.576));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.4}, ParamCoord{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.192));
}

void mock_parameterSpace_nurbs(const std::shared_ptr<NiceMock<MockParameterSpace1>> &parameter_space) {
  set_get_basis_function_nurbs(parameter_space);
  set_basis_function_derivative1(parameter_space);
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(_))
      .WillByDefault(Return(std::array<int, 2>{0, 0}));
  ON_CALL(*parameter_space, GetDegree(_))
      .WillByDefault(Return(Degree{2}));
}

void mock_weights(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace1>> &w_physical_space) {
  ON_CALL(*w_physical_space, GetWeight(_))
      .WillByDefault(Return(1));
  ON_CALL(*w_physical_space, GetWeight(std::array<int, 2>{1, 2}))
      .WillByDefault(Return(2));
}

void mock_homogenous(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace1>> &w_physical_space) {
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{0, 0}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 0.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{1, 0}))
      .WillByDefault(Return(baf::ControlPoint({1.0, 0.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{2, 0}))
      .WillByDefault(Return(baf::ControlPoint({3.0, 0.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{0, 1}))
      .WillByDefault(Return(baf::ControlPoint({-1.0, 0.5})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{1, 1}))
      .WillByDefault(Return(baf::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{2, 1}))
      .WillByDefault(Return(baf::ControlPoint({4.0, 1.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{0, 2}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{1, 2}))
      .WillByDefault(Return(baf::ControlPoint({5.0, 7.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{2, 2}))
      .WillByDefault(Return(baf::ControlPoint({5.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{0, 0}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 0.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{1, 0}))
      .WillByDefault(Return(baf::ControlPoint({1.0, 0.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{2, 0}))
      .WillByDefault(Return(baf::ControlPoint({3.0, 0.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{0, 1}))
      .WillByDefault(Return(baf::ControlPoint({-1.0, 0.5})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{1, 1}))
      .WillByDefault(Return(baf::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{2, 1}))
      .WillByDefault(Return(baf::ControlPoint({4.0, 1.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{0, 2}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{1, 2}))
      .WillByDefault(Return(baf::ControlPoint({2.5, 3.5})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{2, 2}))
      .WillByDefault(Return(baf::ControlPoint({5.0, 2.0})));
}

void mock_weightedPhysicalSpace(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace1>>
                                &w_physical_space) {
  mock_weights(w_physical_space);
  mock_homogenous(w_physical_space);
}

/* 2-dimensional nurbs spline with following properties :
 * KnotVector = {{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}
 * ControlPoints = {{0, 0}, {1, 0}, {3, 0}, {-1, 0.5}, {2, 2}, {4, 1}, {0, 2}, {2.5, 3.5}, {5, 2}}
 * Weights = {1, 1, 1, 1, 1, 1, 1, 2, 1}
*/

class A2DNurbs : public Test {
 public:
  A2DNurbs() :
      parameter_space(std::make_shared<NiceMock<MockParameterSpace1>>()),
      w_physical_space(std::make_shared<NiceMock<MockWeightedPhysicalSpace1>>()){
    /*
    std::array<baf::KnotVector, 2> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}),
         baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})};
    std::array<Degree, 2> degree = {Degree{2}, Degree{2}};
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
    std::array<int, 2> number_of_points = {3, 3};
    std::array<std::shared_ptr<baf::KnotVector>, 2>
        knot_vector_ptr =
        {std::make_shared<baf::KnotVector>(knot_vector[0]), std::make_shared<baf::KnotVector>(knot_vector[1])};

    parameter_space = std::make_shared<spl::ParameterSpace<2>>(spl::ParameterSpace<2>(knot_vector_ptr, degree));
    w_physical_space = std::make_shared<spl::WeightedPhysicalSpace<2>>(spl::WeightedPhysicalSpace<2>(control_points, weights, number_of_points));
    */
    spl::NURBSGenerator<2> nurbs_generator(w_physical_space, parameter_space);
    nurbs_ = std::make_unique<spl::NURBS<2>>(nurbs_generator);
  }

 protected:
  std::unique_ptr<spl::NURBS<2>> nurbs_;
  std::shared_ptr<NiceMock<MockParameterSpace1>> parameter_space;
  std::shared_ptr<NiceMock<MockWeightedPhysicalSpace1>> w_physical_space;
};

TEST_F(A2DNurbs, Returns1_6For0_4And0_6AndDim0) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{0.4}, ParamCoord{0.6}}, {0})[0], DoubleNear(1.62074, 0.00001));
}

TEST_F(A2DNurbs, Returns1_9For0_4And0_6AndDim1) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{0.4}, ParamCoord{0.6}}, {1})[0], DoubleNear(1.88267, 0.00001));
}

TEST_F(A2DNurbs, Returns2_5For0_5And1_0AndDim0) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{0.5}, ParamCoord{1.0}}, {0})[0],
              DoubleNear(2.5, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(A2DNurbs, Returns3_0For0_5And1_0AndDim1) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{0.5}, ParamCoord{1.0}}, {1})[0],
              DoubleNear(3.0, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(A2DNurbs, Returns4_2For0_9And1_0AndDim0) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{0.9}, ParamCoord{1.0}}, {0})[0], DoubleNear(4.19492, 0.00001));
}

TEST_F(A2DNurbs, Returns2_5For0_9And1_0AndDim1) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{0.9}, ParamCoord{1.0}}, {1})[0], DoubleNear(2.45763, 0.00001));
}

TEST_F(A2DNurbs, EvaluatesMultipleValues) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{1.0}, ParamCoord{0.0}}, {0, 1})[0], DoubleEq(3.0));
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{1.0}, ParamCoord{0.0}}, {0, 1})[1], DoubleEq(0.0));
}

TEST_F(A2DNurbs, Returns10_0For0_0And1_0ForDerivative1And0AndDim0) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{1.0}}, {0}, {1, 0})[0], DoubleEq(10.0));
}

TEST_F(A2DNurbs, Returns6_0For0_0And1_0ForDerivative1And0AndDim1) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{1.0}}, {1}, {1, 0})[0], DoubleEq(6.0));
}

TEST_F(A2DNurbs, Returns2_0For0_0And1_0ForDerivative0And1AndDim0) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{1.0}}, {0}, {0, 1})[0], DoubleEq(2.0));
}

TEST_F(A2DNurbs, Returns3_0For0_0And1_0ForDerivative0And1AndDim1) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{1.0}}, {1}, {0, 1})[0], DoubleEq(3.0));
}

TEST_F(A2DNurbs, Returns4_2For0_4And0_6ForDerivative1And0AndDim0) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.4}, ParamCoord{0.6}}, {0}, {1, 0})[0],
              DoubleNear(4.15298, 0.000001));
}

TEST_F(A2DNurbs, Returns0_8For0_4And0_6ForDerivative1And0AndDim1) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.4}, ParamCoord{0.6}}, {1}, {1, 0})[0],
              DoubleNear(0.792032, 0.000001));
}

TEST_F(A2DNurbs, Returns1_4For0_4And0_6ForDerivative0And1AndDim0) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.4}, ParamCoord{0.6}}, {0}, {0, 1})[0],
              DoubleNear(1.40046, 0.00001));
}

TEST_F(A2DNurbs, Returns3_1For0_4And0_6ForDerivative0And1AndDim1) { // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.4}, ParamCoord{0.6}}, {1}, {0, 1})[0],
              DoubleNear(3.13402, 0.00001));
}

class A2DNurbsWithAllWeights1 : public Test {
 public:
  A2DNurbsWithAllWeights1() {
    std::array<baf::KnotVector, 2> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{2}, ParamCoord{2}, ParamCoord{2}}),
         baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{2}, ParamCoord{2}, ParamCoord{2}})};
    std::array<Degree, 2> degree = {Degree{2}, Degree{2}};
    std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 2.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.5})),
        baf::ControlPoint(std::vector<double>({3.0, 4.0})),
        baf::ControlPoint(std::vector<double>({5.0, 3.0})),
        baf::ControlPoint(std::vector<double>({1.0, 4.0})),
        baf::ControlPoint(std::vector<double>({3.5, 5.5})),
        baf::ControlPoint(std::vector<double>({6.0, 4.0}))
    };
    std::array<std::shared_ptr<baf::KnotVector>, 2>
        knot_vector_ptr =
        {std::make_shared<baf::KnotVector>(knot_vector[0]), std::make_shared<baf::KnotVector>(knot_vector[1])};
    nurbs_ = std::make_unique<spl::NURBS<2>>(knot_vector_ptr, degree, control_points, weights);
    bspline_ = std::make_unique<spl::BSpline<2>>(knot_vector_ptr, degree, control_points);
  }

 protected:
  std::unique_ptr<spl::NURBS<2>> nurbs_;
  std::unique_ptr<spl::BSpline<2>> bspline_;
};

TEST_F(A2DNurbsWithAllWeights1, ReturnsSameDerivativeAs2DBSplineFor0_5And0_5AndDerivatives1And1) { // NOLINT
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}}, {0}, {1, 1})[0],
              DoubleEq(bspline_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}}, {0}, {1, 1})[0]));
}

TEST_F(A2DNurbsWithAllWeights1, ReturnsSameDerivativeAs2DBSplineFor0_0And0_7AndDerivatives1And1) { // NOLINT
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.7}}, {0}, {1, 1})[0],
              DoubleEq(bspline_->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.7}}, {0}, {1, 1})[0]));
}

TEST_F(A2DNurbsWithAllWeights1, ReturnsSameDerivativeAs2DBSplineFor0_0And0_7AndDerivatives2And1) { // NOLINT
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.7}}, {0}, {2, 1})[0],
              DoubleNear(bspline_->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.7}}, {0}, {2, 1})[0], 0.000001));
}
