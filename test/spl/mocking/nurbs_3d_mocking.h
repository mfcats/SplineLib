/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef TEST_SPL_MOCKING_NURBS_3D_MOCKING_H_
#define TEST_SPL_MOCKING_NURBS_3D_MOCKING_H_

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

class MockParameterSpace3d : public spl::ParameterSpace<3> {
 public:
  MOCK_CONST_METHOD1(GetDegree, Degree(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 3>, std::array<ParametricCoordinate, 3>));
  MOCK_CONST_METHOD3(GetBasisFunctionDerivatives,
                     double(std::array<int, 3>, std::array<ParametricCoordinate, 3>, std::array<int, 3>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 3>(std::array<ParametricCoordinate, 3>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParametricCoordinate, 3>));
};

class MockWeightedPhysicalSpace3d : public spl::WeightedPhysicalSpace<3> {
 public:
  MOCK_CONST_METHOD1(GetWeight, Weight(std::array<int, 3> const &));
  MOCK_CONST_METHOD1(GetHomogenousControlPoint, spl::ControlPoint(std::array<int, 3>));
  MOCK_CONST_METHOD1(GetControlPoint, spl::ControlPoint(std::array<int, 3> const &));
  MOCK_CONST_METHOD0(GetDimensionality, int());
};

class MockPhysicalSpace3d : public spl::PhysicalSpace<3> {
 public:
  MOCK_CONST_METHOD1(GetControlPoint, spl::ControlPoint(std::array<int, 3> const &));
};

void mock_weights3d(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace3d>> &w_physical_space) {
  ON_CALL(*w_physical_space, GetWeight(_))
      .WillByDefault(Return(Weight{1.0}));
}

void mock_homogenous3d(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace3d>> &w_physical_space) {
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{0, 0, 0}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 0.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{1, 0, 0}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 0.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{2, 0, 0}))
      .WillByDefault(Return(spl::ControlPoint({3.0, 0.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{0, 1, 0}))
      .WillByDefault(Return(spl::ControlPoint({-1.0, 0.5})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(spl::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{2, 1, 0}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 1.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{0, 0, 1}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{1, 0, 1}))
      .WillByDefault(Return(spl::ControlPoint({-1.0, 0.5})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{2, 0, 1}))
      .WillByDefault(Return(spl::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{0, 1, 1}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 1.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{2, 1, 1}))
      .WillByDefault(Return(spl::ControlPoint({5.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{0, 0, 0}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 0.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{1, 0, 0}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 0.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{2, 0, 0}))
      .WillByDefault(Return(spl::ControlPoint({3.0, 0.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{0, 1, 0}))
      .WillByDefault(Return(spl::ControlPoint({-1.0, 0.5})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(spl::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{2, 1, 0}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 1.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{0, 0, 1}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{1, 0, 1}))
      .WillByDefault(Return(spl::ControlPoint({-1.0, 0.5})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{2, 0, 1}))
      .WillByDefault(Return(spl::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{0, 1, 1}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 1.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{2, 1, 1}))
      .WillByDefault(Return(spl::ControlPoint({5.0, 2.0})));
}

void mock_weightedPhysicalSpace3d(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace3d>> &w_physical_space) {
  mock_weights3d(w_physical_space);
  mock_homogenous3d(w_physical_space);
  ON_CALL(*w_physical_space, GetDimensionality()).WillByDefault(Return(2));
}

void mock_physicalSpace3d(const std::shared_ptr<NiceMock<MockPhysicalSpace3d>> &physical_space) {
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{0, 0, 0}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{1, 0, 0}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{2, 0, 0}))
      .WillByDefault(Return(spl::ControlPoint({3.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{0, 1, 0}))
      .WillByDefault(Return(spl::ControlPoint({-1.0, 0.5})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(spl::ControlPoint({2.0, 2.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{2, 1, 0}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 1.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{0, 0, 1}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 2.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{1, 0, 1}))
      .WillByDefault(Return(spl::ControlPoint({-1.0, 0.5})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{2, 0, 1}))
      .WillByDefault(Return(spl::ControlPoint({2.0, 2.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{0, 1, 1}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 1.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 2.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{2, 1, 1}))
      .WillByDefault(Return(spl::ControlPoint({5.0, 2.0})));
}

void set_get_basis_function_nurbs3d(const std::shared_ptr<NiceMock<MockParameterSpace3d>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctions(_,
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.0625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 0, 0},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 1, 0},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 0, 1},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 1, 1},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{0, 0, 0},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.8},
                                                                                  ParametricCoordinate{0.1}}))
      .WillByDefault(Return(0.045));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 0, 0},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.8},
                                                                                  ParametricCoordinate{0.1}}))
      .WillByDefault(Return(0.09));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{2, 0, 0},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.8},
                                                                                  ParametricCoordinate{0.1}}))
      .WillByDefault(Return(0.045));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{0, 1, 0},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.8},
                                                                                  ParametricCoordinate{0.1}}))
      .WillByDefault(Return(0.18));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 1, 0},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.8},
                                                                                  ParametricCoordinate{0.1}}))
      .WillByDefault(Return(0.36));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{2, 1, 0},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.8},
                                                                                  ParametricCoordinate{0.1}}))
      .WillByDefault(Return(0.18));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{0, 0, 1},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.8},
                                                                                  ParametricCoordinate{0.1}}))
      .WillByDefault(Return(0.005));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 0, 1},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.8},
                                                                                  ParametricCoordinate{0.1}}))
      .WillByDefault(Return(0.01));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{2, 0, 1},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.8},
                                                                                  ParametricCoordinate{0.1}}))
      .WillByDefault(Return(0.005));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{0, 1, 1},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.8},
                                                                                  ParametricCoordinate{0.1}}))
      .WillByDefault(Return(0.02));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 1, 1},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.8},
                                                                                  ParametricCoordinate{0.1}}))
      .WillByDefault(Return(0.04));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{2, 1, 1},
                                              std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.8},
                                                                                  ParametricCoordinate{0.1}}))
      .WillByDefault(Return(0.02));
}

void set_basis_function_derivative3d(const std::shared_ptr<NiceMock<MockParameterSpace3d>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(_, _, _)).WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 0, 0},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 0, 0},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(-0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 1, 0},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(-0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 1, 0},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 0, 1},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 0, 1},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(-0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 1, 1},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(-0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 1, 1},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 0, 0},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.8},
                                                                                            ParametricCoordinate{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(-1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 0, 0},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.8},
                                                                                            ParametricCoordinate{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 1, 0},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.8},
                                                                                            ParametricCoordinate{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 1, 0},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.8},
                                                                                            ParametricCoordinate{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(-1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 0, 1},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.8},
                                                                                            ParametricCoordinate{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 0, 1},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.8},
                                                                                            ParametricCoordinate{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(-1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 1, 1},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.8},
                                                                                            ParametricCoordinate{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(-1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 1, 1},
                                                        std::array<ParametricCoordinate, 3>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.8},
                                                                                            ParametricCoordinate{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(1.0));
}

void mock_parameterSpace_nurbs3d(const std::shared_ptr<NiceMock<MockParameterSpace3d>> &parameter_space) {
  set_get_basis_function_nurbs3d(parameter_space);
  set_basis_function_derivative3d(parameter_space);
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(_))
      .WillByDefault(Return(std::array<int, 3>{0, 0, 0}));
  ON_CALL(*parameter_space, GetDegree(_))
      .WillByDefault(Return(Degree{1}));
  ON_CALL(*parameter_space, GetDegree(0))
      .WillByDefault(Return(Degree{2}));
}

#endif  // TEST_SPL_MOCKING_NURBS_3D_MOCKING_H_
