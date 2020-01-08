/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef TEST_SPL_MOCKING_NURBS_1D_MOCKING_H_
#define TEST_SPL_MOCKING_NURBS_1D_MOCKING_H_

#include "gmock/gmock.h"

#include "src/spl/nurbs.h"
#include "src/util/numeric_settings.h"

using ::testing::NiceMock;
using testing::Return;
using ::testing::Throw;
using ::testing::_;

using namespace splinelib::src;

class MockParameterSpace14111 : public spl::ParameterSpace<1> {
 public:
  MOCK_CONST_METHOD1(GetDegree, Degree(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 1>, std::array<ParametricCoordinate, 1>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 1>(std::array<ParametricCoordinate, 1>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParametricCoordinate, 1>));
};

class MockWeightedPhysicalSpace14111 : public spl::WeightedPhysicalSpace<1> {
 public:
  MOCK_CONST_METHOD1(GetWeight, Weight(std::array<int, 1> const &));
  MOCK_CONST_METHOD1(GetHomogenousControlPoint, spl::ControlPoint(std::array<int, 1>));
  MOCK_CONST_METHOD0(GetDimensionality, int());
};

void mock_weights(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace14111>> &w_physical_space) {
  ON_CALL(*w_physical_space, GetWeight(_))
      .WillByDefault(Return(Weight{1.0}));
  ON_CALL(*w_physical_space, GetWeight(std::array<int, 1>{1}))
      .WillByDefault(Return(Weight{4.0}));
}

void mock_homogenous(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace14111>> &w_physical_space) {
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 1>{0}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 0.0, 1.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 1>{1}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 4.0, 4.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 1>{2}))
      .WillByDefault(Return(spl::ControlPoint({3.0, 2.0, 1.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 1>{3}))
      .WillByDefault(Return(spl::ControlPoint({4.0, 1.0, 1.0})));
}

void mock_weightedPhysicalSpace(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace14111>>
                                &w_physical_space) {
  mock_weights(w_physical_space);
  mock_homogenous(w_physical_space);
  ON_CALL(*w_physical_space, GetDimensionality()).WillByDefault(Return(2));
}

void set_get_basis_function_nurbs(const std::shared_ptr<NiceMock<MockParameterSpace14111>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctions(_, std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0}}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{3},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0}}))
      .WillByDefault(Return(0));
}

void mock_parameterSpace_nurbs(const std::shared_ptr<NiceMock<MockParameterSpace14111>> &parameter_space) {
  set_get_basis_function_nurbs(parameter_space);
  ON_CALL(*parameter_space,
          GetArrayOfFirstNonZeroBasisFunctions(std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0}}))
      .WillByDefault(Return(std::array<int, 1>{1}));
  ON_CALL(*parameter_space, GetDegree(0))
      .WillByDefault(Return(Degree{2}));
}

class MockParameterSpace1009 : public spl::ParameterSpace<1> {
 public:
  MOCK_CONST_METHOD1(GetDegree, Degree(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 1>, std::array<ParametricCoordinate, 1>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 1>(std::array<ParametricCoordinate, 1>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParametricCoordinate, 1>));
};

class MockWeightedPhysicalSpace1009 : public spl::WeightedPhysicalSpace<1> {
 public:
  MOCK_CONST_METHOD1(GetWeight, Weight(std::array<int, 1> const &));
  MOCK_CONST_METHOD1(GetHomogenousControlPoint, spl::ControlPoint(std::array<int, 1>));
  MOCK_CONST_METHOD0(GetDimensionality, int());
};

void mock_weights(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace1009>> &w_physical_space) {
  ON_CALL(*w_physical_space, GetWeight(std::array<int, 1>{0}))
      .WillByDefault(Return(Weight{1.0}));
  ON_CALL(*w_physical_space, GetWeight(std::array<int, 1>{1}))
      .WillByDefault(Return(Weight{0.9}));
  ON_CALL(*w_physical_space, GetWeight(std::array<int, 1>{2}))
      .WillByDefault(Return(Weight{0.7}));
  ON_CALL(*w_physical_space, GetWeight(std::array<int, 1>{3}))
      .WillByDefault(Return(Weight{0.5}));
  ON_CALL(*w_physical_space, GetWeight(std::array<int, 1>{4}))
      .WillByDefault(Return(Weight{0.8}));
  ON_CALL(*w_physical_space, GetWeight(std::array<int, 1>{5}))
      .WillByDefault(Return(Weight{1.2}));
  ON_CALL(*w_physical_space, GetWeight(std::array<int, 1>{6}))
      .WillByDefault(Return(Weight{2.0}));
}

void mock_homogenous(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace1009>> &w_physical_space) {
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 1>{0}))
      .WillByDefault(Return(spl::ControlPoint({0.5, 3.0, 1.0, 1.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 1>{1}))
      .WillByDefault(Return(spl::ControlPoint({1.35, 4.95, 3.6, 0.9})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 1>{2}))
      .WillByDefault(Return(spl::ControlPoint({3.15, 3.85, 0.07, 0.7})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 1>{3}))
      .WillByDefault(Return(spl::ControlPoint({1.5, 0.75, 1.0, 0.5})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 1>{4}))
      .WillByDefault(Return(spl::ControlPoint({6.0, 1.2, 2.8, 0.8})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 1>{5}))
      .WillByDefault(Return(spl::ControlPoint({7.2, 4.8, 6.36, 1.2})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 1>{6}))
      .WillByDefault(Return(spl::ControlPoint({17.0, 9.0, 0.0, 2.0})));
}

void mock_weightedPhysicalSpace(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace1009>>
                                &w_physical_space) {
  mock_weights(w_physical_space);
  mock_homogenous(w_physical_space);
  ON_CALL(*w_physical_space, GetDimensionality()).WillByDefault(Return(3));
}

void set_throw_method(const std::shared_ptr<NiceMock<MockParameterSpace1009>> &parameter_space) {
  ON_CALL(*parameter_space,
          ThrowIfParametricCoordinateOutsideKnotVectorRange(std::array<ParametricCoordinate, 1>{
              ParametricCoordinate{1.2}}))
      .WillByDefault(Throw(std::range_error("Out of knotvector range")));
  ON_CALL(*parameter_space,
          ThrowIfParametricCoordinateOutsideKnotVectorRange(std::array<ParametricCoordinate, 1>{
              ParametricCoordinate{-0.1}}))
      .WillByDefault(Throw(std::range_error("Out of knotvector range")));
}

void set_get_basis_function_nurbs(const std::shared_ptr<NiceMock<MockParameterSpace1009>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctions(_, std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.0}}))
      .WillByDefault(Return(0));
  ON_CALL(*parameter_space, GetBasisFunctions(_, std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.25}}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctions(_, std::array<ParametricCoordinate, 1>{ParametricCoordinate{1}}))
      .WillByDefault(Return(0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{0},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.0}}))
      .WillByDefault(Return(1));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{3},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.25}}))
      .WillByDefault(Return(0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{1},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0 / 3.0}}))
      .WillByDefault(Return(2.0 / 9.0));
}

void set_get_basis_function_nurbs_2(const std::shared_ptr<NiceMock<MockParameterSpace1009>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{2},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0 / 3.0}}))
      .WillByDefault(Return(65.0 / 90.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{3},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0 / 3.0}}))
      .WillByDefault(Return(5.0 / 90.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{1},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0 / 3.0}}))
      .WillByDefault(Return(2.0 / 9.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{2},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0 / 3.0}}))
      .WillByDefault(Return(65.0 / 90.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{3},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0 / 3.0}}))
      .WillByDefault(Return(5.0 / 90.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{6},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{1}}))
      .WillByDefault(Return(1));
}

void mock_parameterSpace_nurbs(const std::shared_ptr<NiceMock<MockParameterSpace1009>> &parameter_space) {
  set_get_basis_function_nurbs(parameter_space);
  set_get_basis_function_nurbs_2(parameter_space);
  set_throw_method(parameter_space);
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(
      std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.0}}))
      .WillByDefault(Return(std::array<int, 1>{0}));
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(
      std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.25}}))
      .WillByDefault(Return(std::array<int, 1>{1}));
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(
      std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0 / 3.0}}))
      .WillByDefault(Return(std::array<int, 1>{1}));
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(
      std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0}}))
      .WillByDefault(Return(std::array<int, 1>{4}));
  ON_CALL(*parameter_space, GetDegree(0))
      .WillByDefault(Return(Degree{2}));
}

class MockParameterSpace112 : public spl::ParameterSpace<1> {
 public:
  MOCK_CONST_METHOD1(GetDegree, Degree(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 1>, std::array<ParametricCoordinate, 1>));
  MOCK_CONST_METHOD3(GetBasisFunctionDerivatives,
                     double(std::array<int, 1>, std::array<ParametricCoordinate, 1>, std::array<int, 1>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 1>(std::array<ParametricCoordinate, 1>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParametricCoordinate, 1>));
};

class MockWeightedPhysicalSpace112 : public spl::WeightedPhysicalSpace<1> {
 public:
  MOCK_CONST_METHOD1(GetControlPoint, spl::ControlPoint(std::array<int, 1> const &));
  MOCK_CONST_METHOD1(GetWeight, Weight(std::array<int, 1> const &));
  MOCK_CONST_METHOD0(GetDimensionality, int());
};

void mock_weights(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace112>> &w_physical_space) {
  ON_CALL(*w_physical_space, GetWeight(_))
      .WillByDefault(Return(Weight{1.0}));
  ON_CALL(*w_physical_space, GetWeight(std::array<int, 1>{2}))
      .WillByDefault(Return(Weight{2.0}));
}

void mock_weightedPhysicalSpace(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace112>>
                                &w_physical_space) {
  mock_weights(w_physical_space);
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 1>{0}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 0.0, 1.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 1>{1}))
      .WillByDefault(Return(spl::ControlPoint({1.0, 1.0, 1.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 1>{2}))
      .WillByDefault(Return(spl::ControlPoint({0.0, 1.0, 2.0})));
  ON_CALL(*w_physical_space, GetDimensionality()).WillByDefault(Return(2));
}

void set_get_basis_function_nurbs(const std::shared_ptr<NiceMock<MockParameterSpace112>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctions(_, std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.0}}))
      .WillByDefault(Return(0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{0},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.0}}))
      .WillByDefault(Return(1));
  ON_CALL(*parameter_space, GetBasisFunctions(_,
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.25));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{1},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.5}}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctions(_,
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{1}}))
      .WillByDefault(Return(0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{2},
                                              std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0}}))
      .WillByDefault(Return(1));
}

void set_basis_function_derivative1_nurbs(const std::shared_ptr<NiceMock<MockParameterSpace112>> &parameter_space) {
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(_, _, _)).WillByDefault(Return(0.0));
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{0},
                                      std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.0}},
                                      std::array<int, 1>{1})).WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{1},
                                      std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.0}},
                                      std::array<int, 1>{1})).WillByDefault(Return(2.0));
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{0},
                                      std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.5}},
                                      std::array<int, 1>{1})).WillByDefault(Return(-1.0));
}

void set_basis_function_derivative1_nurbs_2(const std::shared_ptr<NiceMock<MockParameterSpace112>> &parameter_space) {
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{2},
                                      std::array<ParametricCoordinate, 1>{ParametricCoordinate{0.5}},
                                      std::array<int, 1>{1})).WillByDefault(Return(1.0));
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{1},
                                      std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0}},
                                      std::array<int, 1>{1})).WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{2},
                                      std::array<ParametricCoordinate, 1>{ParametricCoordinate{1.0}},
                                      std::array<int, 1>{1})).WillByDefault(Return(2.0));
}

void set_basis_function_derivative2_nurbs(const std::shared_ptr<NiceMock<MockParameterSpace112>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1>{0}, _,
                                                        std::array<int, 1>{2})).WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1>{1}, _,
                                                        std::array<int, 1>{2})).WillByDefault(Return(-4.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1>{2}, _,
                                                        std::array<int, 1>{2})).WillByDefault(Return(2.0));
}

void mock_parameterSpace_nurbs(const std::shared_ptr<NiceMock<MockParameterSpace112>> &parameter_space) {
  set_get_basis_function_nurbs(parameter_space);
  set_basis_function_derivative1_nurbs(parameter_space);
  set_basis_function_derivative1_nurbs_2(parameter_space);
  set_basis_function_derivative2_nurbs(parameter_space);
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(_))
      .WillByDefault(Return(std::array<int, 1>{0}));
  ON_CALL(*parameter_space, GetDegree(0))
      .WillByDefault(Return(Degree{2}));
}

#endif  // TEST_SPL_MOCKING_NURBS_1D_MOCKING_H_
