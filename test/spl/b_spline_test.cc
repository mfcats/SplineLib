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

using testing::Test;
using testing::DoubleEq;
using ::testing::Return;
using ::testing::Throw;
using ::testing::NiceMock;

class MockParameterSpace : public spl::ParameterSpace<1> {
 public:
  MOCK_CONST_METHOD1(GetDegree, Degree(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 1>, std::array<ParamCoord, 1>));
  MOCK_CONST_METHOD3(GetBasisFunctionDerivatives,
                     double(std::array<int, 1>, std::array<ParamCoord, 1>, std::array<int, 1>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 1>(std::array<ParamCoord, 1>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParamCoord, 1>));
};

class MockPhysicalSpace : public spl::PhysicalSpace<1> {
 public:
  MOCK_CONST_METHOD1(GetControlPoint, baf::ControlPoint(std::array<int, 1>));
};

void mock_physicalSpace(const std::shared_ptr<NiceMock<MockPhysicalSpace>> &physical_space) {
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{0}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{1}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 1.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{2}))
      .WillByDefault(Return(baf::ControlPoint({1.0, 1.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{3}))
      .WillByDefault(Return(baf::ControlPoint({1.5, 1.5})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{4}))
      .WillByDefault(Return(baf::ControlPoint({2.0, 1.3})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{5}))
      .WillByDefault(Return(baf::ControlPoint({3.0, 2.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{6}))
      .WillByDefault(Return(baf::ControlPoint({4.0, 1.5})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 1>{7}))
      .WillByDefault(Return(baf::ControlPoint({4.0, 0.0})));
}

void set_throw_method(const std::shared_ptr<NiceMock<MockParameterSpace>> &parameter_space) {
  ON_CALL(*parameter_space,
          ThrowIfParametricCoordinateOutsideKnotVectorRange(std::array<ParamCoord, 1>{ParamCoord{-1.0}}))
      .WillByDefault(Throw(std::range_error("Out of knotvector range")));
  ON_CALL(*parameter_space,
          ThrowIfParametricCoordinateOutsideKnotVectorRange(std::array<ParamCoord, 1>{ParamCoord{6.0}}))
      .WillByDefault(Throw(std::range_error("Out of knotvector range")));
}

void set_get_basis_function(const std::shared_ptr<NiceMock<MockParameterSpace>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{0}, std::array<ParamCoord, 1>{ParamCoord{0.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{2}, std::array<ParamCoord, 1>{ParamCoord{2.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{3}, std::array<ParamCoord, 1>{ParamCoord{2.5}}))
      .WillByDefault(Return(0.75));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{4}, std::array<ParamCoord, 1>{ParamCoord{2.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{5}, std::array<ParamCoord, 1>{ParamCoord{5.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{6}, std::array<ParamCoord, 1>{ParamCoord{5.0}}))
      .WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 1>{7}, std::array<ParamCoord, 1>{ParamCoord{5.0}}))
      .WillByDefault(Return(1.0));
}

void set_basis_function_derivative1(const std::shared_ptr<NiceMock<MockParameterSpace>> &parameter_space) {
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{0},
                                      std::array<ParamCoord, 1>{ParamCoord{0.0}},
                                      std::array<int, 1>{1})).WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{1},
                                      std::array<ParamCoord, 1>{ParamCoord{0.0}},
                                      std::array<int, 1>{1})).WillByDefault(Return(2.0));
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{2},
                                      std::array<ParamCoord, 1>{ParamCoord{0.0}},
                                      std::array<int, 1>{1})).WillByDefault(Return(0.0));
  ON_CALL(*parameter_space,
          GetBasisFunctionDerivatives(std::array<int, 1>{5},
                                      std::array<ParamCoord, 1>{ParamCoord{5.0}},
                                      std::array<int, 1>{1})).WillByDefault(Return(0.0));
}

void set_basis_function_derivative2(const std::shared_ptr<NiceMock<MockParameterSpace>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1>{6},
                                                        std::array<ParamCoord, 1>{ParamCoord{5.0}},
                                                        std::array<int, 1>{1})).WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1>{7},
                                                        std::array<ParamCoord, 1>{ParamCoord{5.0}},
                                                        std::array<int, 1>{1})).WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1>{2},
                                                        std::array<ParamCoord, 1>{ParamCoord{2.25}},
                                                        std::array<int, 1>{1})).WillByDefault(Return(-0.75));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1>{3},
      std::array<ParamCoord, 1>{ParamCoord{2.25}},
                                      std::array<int, 1>{1})).WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 1>{4},
                                      std::array<ParamCoord, 1>{ParamCoord{2.25}},
                                      std::array<int, 1>{1})).WillByDefault(Return(0.25));
}

void mock_parameterSpace(const std::shared_ptr<NiceMock<MockParameterSpace>> &parameter_space) {
  set_throw_method(parameter_space);
  set_get_basis_function(parameter_space);
  set_basis_function_derivative1(parameter_space);
  set_basis_function_derivative2(parameter_space);
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(std::array<ParamCoord, 1>{ParamCoord{0.0}}))
      .WillByDefault(Return(std::array<int, 1>{0}));
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(std::array<ParamCoord, 1>{ParamCoord{2.25}}))
      .WillByDefault(Return(std::array<int, 1>{2}));
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(std::array<ParamCoord, 1>{ParamCoord{2.5}}))
      .WillByDefault(Return(std::array<int, 1>{2}));
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(std::array<ParamCoord, 1>{ParamCoord{5.0}}))
      .WillByDefault(Return(std::array<int, 1>{5}));
  ON_CALL(*parameter_space, GetDegree(0))
      .WillByDefault(Return(Degree{2}));
}

// U = {0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}
class ABSpline : public Test {
 public:
  ABSpline() :
      parameter_space(std::make_shared<NiceMock<MockParameterSpace>>()),
      physical_space(std::make_shared<NiceMock<MockPhysicalSpace>>()) {
    spl::BSplineGenerator<1> b_spline_generator(physical_space, parameter_space);
    b_spline = std::make_unique<spl::BSpline<1>>(b_spline_generator);
  }

 protected:
  std::unique_ptr<spl::BSpline<1>> b_spline;
  std::shared_ptr<NiceMock<MockParameterSpace>> parameter_space;
  std::shared_ptr<NiceMock<MockPhysicalSpace>> physical_space;
};

TEST_F(ABSpline, Returns0_0For0AndDim0) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{0.0}}, {0})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns0_0For0AndDim1) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{0.0}}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns4_0For5AndDim0) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{5.0}}, {0})[0], DoubleEq(4.0));
}

TEST_F(ABSpline, Returns0_0For5AndDim1) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{5.0}}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns1_5For2_5AndDim0) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{2.5}}, {0})[0], DoubleEq(1.5));
}

TEST_F(ABSpline, Returns0_0For0_0Dim0AndDer1) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{0.0}}, {0}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns1_0For0_0Dim1AndDer1) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{0.0}}, {1}, {1})[0], DoubleEq(2.0));
}

TEST_F(ABSpline, Returns12_0For5_0Dim0AndDer1) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{5.0}}, {0}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns0_325For2_25Dim1AndDer1) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->EvaluateDerivative({ParamCoord{2.25}}, {1}, {1})[0], DoubleEq(0.325));
}

TEST_F(ABSpline, CanBeConstructedWithAPhysicalAndAParameterSpace) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{5.0}}, {0})[0], DoubleEq(4.0));
}

TEST_F(ABSpline, ThrowsExceptionForEvaluationAt6_0) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THROW(b_spline->Evaluate({ParamCoord{6.0}}, {0}), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForEvaluationAtMinus1_0) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THROW(b_spline->Evaluate({ParamCoord{-1.0}}, {0}), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForDerivativeEvaluationAt6_0) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THROW(b_spline->EvaluateDerivative({ParamCoord{6.0}}, {0}, {1}), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForDerivativeEvaluationAtMinus1_0) { // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THROW(b_spline->EvaluateDerivative({ParamCoord{-1.0}}, {0}, {1}), std::runtime_error);
}

class ABSplineWithSplineGenerator : public Test {
 public:
  ABSplineWithSplineGenerator() : degree_{Degree{2}} {
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
    knot_vector_[0] = {std::make_shared<baf::KnotVector>(knot_vector[0])};
    spl::BSplineGenerator<1> b_spline_generator(knot_vector_, degree_, control_points_);
    b_spline = std::make_unique<spl::BSpline<1>>(b_spline_generator);
  }

 protected:
  std::unique_ptr<spl::BSpline<1>> b_spline;
  KnotVectors<1> knot_vector_;
  std::array<Degree, 1> degree_;
  std::vector<baf::ControlPoint> control_points_;
};

TEST_F(ABSplineWithSplineGenerator, Returns0_0For0AndDim0) { // NOLINT
  ASSERT_THAT(b_spline->Evaluate({ParamCoord{0.0}}, {0})[0], DoubleEq(0.0));
}
