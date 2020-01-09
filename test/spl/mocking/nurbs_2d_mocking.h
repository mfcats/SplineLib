/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef TEST_SPL_MOCKING_NURBS_2D_MOCKING_H_
#define TEST_SPL_MOCKING_NURBS_2D_MOCKING_H_

#include <array>
#include <numeric>

#include "gmock/gmock.h"

#include "src/spl/nurbs.h"
#include "src/util/numeric_settings.h"

using testing::Test;
using ::testing::Return;
using ::testing::NiceMock;
using ::testing::_;

using namespace splinelib::src;

class MockParameterSpace1 : public spl::ParameterSpace<2> {
 public:
  MOCK_CONST_METHOD1(GetDegree, Degree(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 2>, std::array<ParametricCoordinate, 2>));
  MOCK_CONST_METHOD3(GetBasisFunctionDerivatives,
                     double(std::array<int, 2>, std::array<ParametricCoordinate, 2>, std::array<int, 2>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 2>(std::array<ParametricCoordinate, 2>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParametricCoordinate, 2>));
};

class MockWeightedPhysicalSpace1 : public spl::WeightedPhysicalSpace<2> {
 public:
  MOCK_CONST_METHOD1(GetWeight, Weight(std::array<int, 2> const &));
  MOCK_CONST_METHOD1(GetHomogenousControlPoint, spl::ControlPoint(std::array<int, 2> const &));
  MOCK_CONST_METHOD1(GetControlPoint, spl::ControlPoint(std::array<int, 2> const &));
  MOCK_CONST_METHOD0(GetDimensionality, int());
};

void set_get_basis_function_nurbs(const std::shared_ptr<NiceMock<MockParameterSpace1>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctions(_,
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{1.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{1.0}}))
      .WillByDefault(Return(0.25));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{1.0}}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{1.0}}))
      .WillByDefault(Return(0.25));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                  ParametricCoordinate{0.6}}))
      .WillByDefault(Return(0.0576));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                  ParametricCoordinate{0.6}}))
      .WillByDefault(Return(0.0768));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                  ParametricCoordinate{0.6}}))
      .WillByDefault(Return(0.0256));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                  ParametricCoordinate{0.6}}))
      .WillByDefault(Return(0.1728));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                  ParametricCoordinate{0.6}}))
      .WillByDefault(Return(0.2304));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                  ParametricCoordinate{0.6}}))
      .WillByDefault(Return(0.0768));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                  ParametricCoordinate{0.6}}))
      .WillByDefault(Return(0.1296));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                  ParametricCoordinate{0.6}}))
      .WillByDefault(Return(0.1728));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                  ParametricCoordinate{0.6}}))
      .WillByDefault(Return(0.0576));
  ON_CALL(*parameter_space, GetBasisFunctions(_,
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.9},
                                                                                  ParametricCoordinate{1.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.9},
                                                                                  ParametricCoordinate{1.0}}))
      .WillByDefault(Return(0.01));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.9},
                                                                                  ParametricCoordinate{1.0}}))
      .WillByDefault(Return(0.18));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.9},
                                                                                  ParametricCoordinate{1.0}}))
      .WillByDefault(Return(0.81));
  ON_CALL(*parameter_space, GetBasisFunctions(_,
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{1.0},
                                                                                  ParametricCoordinate{0.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{1.0},
                                                                                  ParametricCoordinate{0.0}}))
      .WillByDefault(Return(1));
  ON_CALL(*parameter_space, GetBasisFunctions(_,
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                  ParametricCoordinate{1.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                  ParametricCoordinate{1.0}}))
      .WillByDefault(Return(1));
}

void set_basis_function_derivative1(const std::shared_ptr<NiceMock<MockParameterSpace1>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(_, _, _)).WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{1.0}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{1.0}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{1.0}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{1.0}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.192));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.064));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.128));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.576));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.192));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.384));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.432));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.144));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.288));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.288));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.384));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.128));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.144));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.192));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.064));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.432));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.576));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.4},
                                                                                            ParametricCoordinate{0.6}},
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
      .WillByDefault(Return(Weight{1.0}));
  ON_CALL(*w_physical_space, GetWeight(std::array<int, 2>{1, 2}))
      .WillByDefault(Return(Weight{2.0}));
}

void mock_homogenous(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace1>> &w_physical_space) {
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{0, 0}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 0.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{1, 0}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 0.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{2, 0}))
      .WillByDefault(Return(spl::ControlPoint({3.0, 0.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{0, 1}))
      .WillByDefault(Return(spl::ControlPoint({-1.0, 0.5})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{1, 1}))
      .WillByDefault(Return(spl::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{2, 1}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 1.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{0, 2}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{1, 2}))
      .WillByDefault(Return(spl::ControlPoint({5.0, 7.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{2, 2}))
      .WillByDefault(Return(spl::ControlPoint({5.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{0, 0}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 0.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{1, 0}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 0.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{2, 0}))
      .WillByDefault(Return(spl::ControlPoint({3.0, 0.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{0, 1}))
      .WillByDefault(Return(spl::ControlPoint({-1.0, 0.5})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{1, 1}))
      .WillByDefault(Return(spl::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{2, 1}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 1.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{0, 2}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{1, 2}))
      .WillByDefault(Return(spl::ControlPoint({2.5, 3.5})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{2, 2}))
      .WillByDefault(Return(spl::ControlPoint({5.0, 2.0})));
}

void mock_weightedPhysicalSpace(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace1>>
                                &w_physical_space) {
  mock_weights(w_physical_space);
  mock_homogenous(w_physical_space);
  ON_CALL(*w_physical_space, GetDimensionality()).WillByDefault(Return(2));
}
class MockParameterSpace2 : public spl::ParameterSpace<2> {
 public:
  MOCK_CONST_METHOD1(GetDegree, Degree(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 2>, std::array<ParametricCoordinate, 2>));
  MOCK_CONST_METHOD3(GetBasisFunctionDerivatives,
                     double(std::array<int, 2>, std::array<ParametricCoordinate, 2>, std::array<int, 2>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 2>(std::array<ParametricCoordinate, 2>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParametricCoordinate, 2>));
};

class MockPhysicalSpace2 : public spl::PhysicalSpace<2> {
 public:
  MOCK_CONST_METHOD1(GetControlPoint, spl::ControlPoint(std::array<int, 2> const &));
};

class MockWeightedPhysicalSpace2 : public spl::WeightedPhysicalSpace<2> {
 public:
  MOCK_CONST_METHOD1(GetWeight, Weight(std::array<int, 2> const &));
  MOCK_CONST_METHOD1(GetHomogenousControlPoint, spl::ControlPoint(std::array<int, 2> const &));
  MOCK_CONST_METHOD1(GetControlPoint, spl::ControlPoint(std::array<int, 2> const &));
  MOCK_CONST_METHOD0(GetDimensionality, int());
};

void set_get_basis_function_nurbs(const std::shared_ptr<NiceMock<MockParameterSpace2>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctions(_,
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                  ParametricCoordinate{0.7}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.316406));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.210938));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.0351562));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.210938));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.140625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.0234375));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.0351562));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.0234375));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.00390625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                  ParametricCoordinate{0.7}}))
      .WillByDefault(Return(0.4225));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                  ParametricCoordinate{0.7}}))
      .WillByDefault(Return(0.455));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                  ParametricCoordinate{0.7}}))
      .WillByDefault(Return(0.1225));
}

void set_basis_function_derivative1(const std::shared_ptr<NiceMock<MockParameterSpace2>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(0.5625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(-0.375));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(-0.1875));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(-0.375));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(0.25));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(-0.1875));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(0.0625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(0.65));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(-0.65));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(-0.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(-0.3));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(0.3));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(-0.35));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(0.35));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{1, 1}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{2, 1}))
      .WillByDefault(Return(-0.325));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{2, 1}))
      .WillByDefault(Return(0.65));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{2, 1}))
      .WillByDefault(Return(-0.325));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{2, 1}))
      .WillByDefault(Return(0.15));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{2, 1}))
      .WillByDefault(Return(-0.3));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{2, 1}))
      .WillByDefault(Return(0.15));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{2, 1}))
      .WillByDefault(Return(0.175));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{2, 1}))
      .WillByDefault(Return(-0.35));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.7}},
                                                        std::array<int, 2>{2, 1}))
      .WillByDefault(Return(0.175));
}

void mock_parameterSpace_nurbs(const std::shared_ptr<NiceMock<MockParameterSpace2>> &parameter_space) {
  set_get_basis_function_nurbs(parameter_space);
  set_basis_function_derivative1(parameter_space);
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(_))
      .WillByDefault(Return(std::array<int, 2>{0, 0}));
  ON_CALL(*parameter_space, GetDegree(_))
      .WillByDefault(Return(Degree{2}));
}

void mock_weights(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace2>> &w_physical_space) {
  ON_CALL(*w_physical_space, GetWeight(_))
      .WillByDefault(Return(Weight{1.0}));
}

void mock_homogenous(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace2>> &w_physical_space) {
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{0, 0}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{1, 0}))
      .WillByDefault(Return(spl::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{2, 0}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{0, 1}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 2.5})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{1, 1}))
      .WillByDefault(Return(spl::ControlPoint({3.0, 4.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{2, 1}))
      .WillByDefault(Return(spl::ControlPoint({5.0, 3.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{0, 2}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 4.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{1, 2}))
      .WillByDefault(Return(spl::ControlPoint({3.5, 5.5})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 2>{2, 2}))
      .WillByDefault(Return(spl::ControlPoint({6.0, 4.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{0, 0}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{1, 0}))
      .WillByDefault(Return(spl::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{2, 0}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{0, 1}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 2.5})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{1, 1}))
      .WillByDefault(Return(spl::ControlPoint({3.0, 4.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{2, 1}))
      .WillByDefault(Return(spl::ControlPoint({5.0, 3.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{0, 2}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 4.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{1, 2}))
      .WillByDefault(Return(spl::ControlPoint({3.5, 5.5})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 2>{2, 2}))
      .WillByDefault(Return(spl::ControlPoint({6.0, 4.0})));
}

void mock_weightedPhysicalSpace(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace2>>
                                &w_physical_space) {
  mock_weights(w_physical_space);
  mock_homogenous(w_physical_space);
  ON_CALL(*w_physical_space, GetDimensionality()).WillByDefault(Return(2));
}

void mock_physicalSpace(const std::shared_ptr<NiceMock<MockPhysicalSpace2>> &physical_space) {
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{0, 0}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 2.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{1, 0}))
      .WillByDefault(Return(spl::ControlPoint({2.0, 2.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{2, 0}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 2.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{0, 1}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 2.5})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{1, 1}))
      .WillByDefault(Return(spl::ControlPoint({3.0, 4.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{2, 1}))
      .WillByDefault(Return(spl::ControlPoint({5.0, 3.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{0, 2}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 4.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{1, 2}))
      .WillByDefault(Return(spl::ControlPoint({3.5, 5.5})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{2, 2}))
      .WillByDefault(Return(spl::ControlPoint({6.0, 4.0})));
}

#endif  // TEST_SPL_MOCKING_NURBS_2D_MOCKING_H_
