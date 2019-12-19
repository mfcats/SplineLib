/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef TEST_SPL_MOCKING_B_SPLINE_2D_MOCKING_H_
#define TEST_SPL_MOCKING_B_SPLINE_2D_MOCKING_H_

#include <array>
#include <numeric>

#include "gmock/gmock.h"

#include "src/spl/b_spline.h"
#include "src/util/numeric_settings.h"

using testing::Test;
using ::testing::Return;
using ::testing::NiceMock;
using ::testing::_;

using namespace splinelib::src;

class Mock2dParameterSpace : public spl::ParameterSpace<2> {
 public:
  MOCK_CONST_METHOD1(GetDegree, Degree(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 2>, std::array<ParametricCoordinate, 2>));
  MOCK_CONST_METHOD3(GetBasisFunctionDerivatives,
                     double(std::array<int, 2>, std::array<ParametricCoordinate, 2>, std::array<int, 2>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 2>(std::array<ParametricCoordinate, 2>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParametricCoordinate, 2>));
};

class Mock2dPhysicalSpace : public spl::PhysicalSpace<2> {
 public:
  MOCK_CONST_METHOD1(GetControlPoint, spl::ControlPoint(std::array<int, 2>));
};

void mock_2dphysicalSpace(const std::shared_ptr<NiceMock<Mock2dPhysicalSpace>> &physical_space) {
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{0, 0}))
      .WillByDefault(Return(spl::ControlPoint({-1.0, -1.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{1, 0}))
      .WillByDefault(Return(spl::ControlPoint({0.0, -1.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{2, 0}))
      .WillByDefault(Return(spl::ControlPoint({1.0, -1.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{0, 1}))
      .WillByDefault(Return(spl::ControlPoint({-1.0, 0.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{1, 1}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 0.0, 1.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{2, 1}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 0.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{0, 2}))
      .WillByDefault(Return(spl::ControlPoint({-1.0, 1.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{1, 2}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 1.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{2, 2}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 1.0, 0.0})));
}

void set_get_basis_function(const std::shared_ptr<NiceMock<Mock2dParameterSpace>> &parameter_space) {
  // Using wildcards
  ON_CALL(*parameter_space, GetBasisFunctions(_, _))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                  ParametricCoordinate{0.0}}))
      .WillByDefault(Return(1.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                  ParametricCoordinate{0.33333}}))
      .WillByDefault(Return(0.444449));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                  ParametricCoordinate{0.33333}}))
      .WillByDefault(Return(0.444449));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                  ParametricCoordinate{0.33333}}))
      .WillByDefault(Return(0.111109));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.33333},
                                                                                  ParametricCoordinate{0.0}}))
      .WillByDefault(Return(0.444449));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.33333},
                                                                                  ParametricCoordinate{0.0}}))
      .WillByDefault(Return(0.444449));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.33333},
                                                                                  ParametricCoordinate{0.0}}))
      .WillByDefault(Return(0.111109));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.0625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.0625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.25));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.0625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                  ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.0625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                  ParametricCoordinate{0.25}}))
      .WillByDefault(Return(0.0351562));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                  ParametricCoordinate{0.25}}))
      .WillByDefault(Return(0.210938));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 0},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                  ParametricCoordinate{0.25}}))
      .WillByDefault(Return(0.316406));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                  ParametricCoordinate{0.25}}))
      .WillByDefault(Return(0.0234375));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                  ParametricCoordinate{0.25}}))
      .WillByDefault(Return(0.140625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 1},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                  ParametricCoordinate{0.25}}))
      .WillByDefault(Return(0.210938));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                  ParametricCoordinate{0.25}}))
      .WillByDefault(Return(0.00390625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                  ParametricCoordinate{0.25}}))
      .WillByDefault(Return(0.0234375));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 2},
                                              std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                  ParametricCoordinate{0.25}}))
      .WillByDefault(Return(0.0351562));
}

void set_basis_function_derivative1(const std::shared_ptr<NiceMock<Mock2dParameterSpace>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(_, _, _)).WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.0}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.0}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.0}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{0.0}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{
                                                                                                0.33333}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.888898));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{
                                                                                                0.33333}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.888898));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{
                                                                                                0.33333}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.888884));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{
                                                                                                0.33333}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.888884));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{
                                                                                                0.33333}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.222218));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{
                                                                                                0.33333}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.222218));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{
                                                                                                0.33333}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-1.33334));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{
                                                                                                0.33333}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.66668));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.0},
                                                                                            ParametricCoordinate{
                                                                                                0.33333}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.66666));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.25));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.25));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.25));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.5},
                                                                                            ParametricCoordinate{0.5}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.25));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.28125));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.5625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.84375));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.1875));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.375));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.5625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.03125));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.0625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.09375));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.09375));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.5625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.84375));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.0625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.375));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.5625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.03125));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.1875));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.28125));
}

void set_basis_function_derivative2(const std::shared_ptr<NiceMock<Mock2dParameterSpace>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(-1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(3.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(4.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(-6.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(-1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParametricCoordinate, 2>{ParametricCoordinate{0.75},
                                                                                            ParametricCoordinate{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(3.0));
}

void mock_2dparameterSpace(const std::shared_ptr<NiceMock<Mock2dParameterSpace>> &parameter_space) {
  set_basis_function_derivative1(parameter_space);
  set_basis_function_derivative2(parameter_space);
  set_get_basis_function(parameter_space);
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(_))
      .WillByDefault((Return(std::array<int, 2>{0, 0})));
  ON_CALL(*parameter_space, GetDegree(_))
      .WillByDefault(Return(Degree{2}));
}

#endif  // TEST_SPL_MOCKING_B_SPLINE_2D_MOCKING_H_
