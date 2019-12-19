/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include <array>
#include <numeric>

#include "gmock/gmock.h"

#include "src/spl/b_spline.h"

using ::testing::Return;
using ::testing::Throw;
using ::testing::NiceMock;

using namespace splinelib::src;

class MockParameterSpace : public spl::ParameterSpace<1> {
 public:
  MOCK_CONST_METHOD1(GetDegree, Degree(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 1>, std::array<ParametricCoordinate, 1>));
  MOCK_CONST_METHOD3(GetBasisFunctionDerivatives,
                     double(std::array<int, 1>, std::array<ParametricCoordinate, 1>, std::array<int, 1>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 1>(std::array<ParametricCoordinate, 1>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParametricCoordinate, 1>));
};

class MockPhysicalSpace : public spl::PhysicalSpace<1> {
 public:
  MOCK_CONST_METHOD1(GetControlPoint, spl::ControlPoint(std::array<int, 1>));
};

void mock_physicalSpace(const std::shared_ptr<NiceMock<MockPhysicalSpace>> &physical_space) {
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{0}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{1}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 1.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{2}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 1.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{3}))
      .WillByDefault(Return(spl::ControlPoint({1.5, 1.5})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{4}))
      .WillByDefault(Return(spl::ControlPoint({2.0, 1.3})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{5}))
      .WillByDefault(Return(spl::ControlPoint({3.0, 2.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{6}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 1.5})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{7}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 0.0})));
}

void set_throw_method(const std::shared_ptr<NiceMock<MockParameterSpace>> &parameter_space) {
  ON_CALL(*parameter_space,
          ThrowIfParametricCoordinateOutsideKnotVectorRange(std::array<ParametricCoordinate, 1>{
              ParametricCoordinate{-1.0}}))
      .WillByDefault(Throw(std::range_error("Out of knotvector range")));
  ON_CALL(*parameter_space,
          ThrowIfParametricCoordinateOutsideKnotVectorRange(std::array<ParametricCoordinate, 1>{
              ParametricCoordinate{6.0}}))
      .WillByDefault(Throw(std::range_error("Out of knotvector range")));
}

void set_get_basis_function(const std::shared_ptr<NiceMock<MockParameterSpace>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{0},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{2},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{2.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{3},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{2.5}}))
      .WillByDefault(Return(0.75));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{4},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{2.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{5},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{5.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{6},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{5.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{7},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{5.0}}))
      .WillByDefault(Return(1.0));
}

void set_basis_function_derivative1(const std::shared_ptr<NiceMock<MockParameterSpace>> &parameter_space) {
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{0},
                                      std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.0}},
                                      std::array<int, 1>{1})).WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{1},
                                      std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.0}},
                                      std::array<int, 1>{1})).WillByDefault(Return(2.0));
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{2},
                                      std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.0}},
                                      std::array<int, 1>{1})).WillByDefault(Return(0.0));
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{5},
                                      std::array<ParametricCoordinate, 1>{ParametricCoordinate{5.0}},
                                      std::array<int, 1>{1})).WillByDefault(Return(0.0));
}

void set_basis_function_derivative2(const std::shared_ptr<NiceMock<MockParameterSpace>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1>{6},
                                                        std::array<ParametricCoordinate, 1>{ParametricCoordinate{5.0}},
                                                        std::array<int, 1>{1})).WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1>{7},
                                                        std::array<ParametricCoordinate, 1>{ParametricCoordinate{5.0}},
                                                        std::array<int, 1>{1})).WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1>{2},
                                                        std::array<ParametricCoordinate, 1>{ParametricCoordinate{2.25}},
                                                        std::array<int, 1>{1})).WillByDefault(Return(-0.75));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1>{3},
                                                        std::array<ParametricCoordinate, 1>{ParametricCoordinate{2.25}},
                                                        std::array<int, 1>{1})).WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1>{4},
                                                        std::array<ParametricCoordinate, 1>{ParametricCoordinate{2.25}},
                                                        std::array<int, 1>{1})).WillByDefault(Return(0.25));
}

void mock_parameterSpace(const std::shared_ptr<NiceMock<MockParameterSpace>> &parameter_space) {
  set_throw_method(parameter_space);
  set_get_basis_function(parameter_space);
  set_basis_function_derivative1(parameter_space);
  set_basis_function_derivative2(parameter_space);
  ON_CALL(*parameter_space,
          GetArrayOfFirstNonZeroBasisFunctions(std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.0}}))
      .WillByDefault(Return(std::array<int, 1>{0}));
  ON_CALL(*parameter_space,
          GetArrayOfFirstNonZeroBasisFunctions(std::array<ParametricCoordinate, 1>{ParametricCoordinate{2.25}}))
      .WillByDefault(Return(std::array<int, 1>{2}));
  ON_CALL(*parameter_space,
          GetArrayOfFirstNonZeroBasisFunctions(std::array<ParametricCoordinate, 1>{ParametricCoordinate{2.5}}))
      .WillByDefault(Return(std::array<int, 1>{2}));
  ON_CALL(*parameter_space,
          GetArrayOfFirstNonZeroBasisFunctions(std::array<ParametricCoordinate, 1>{ParametricCoordinate{5.0}}))
      .WillByDefault(Return(std::array<int, 1>{5}));
  ON_CALL(*parameter_space, GetDegree(0))
      .WillByDefault(Return(Degree{2}));
}
