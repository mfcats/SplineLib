/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include <array>

#include "gmock/gmock.h"

#include "src/spl/b_spline.h"
#include "test/spl/mocking/b_spline_1d_mocking.h"

using testing::Test;
using testing::DoubleEq;

using namespace splinelib::src;

// U = {0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}
class ABSpline : public Test {
 public:
  ABSpline() :
      parameter_space(std::make_shared<NiceMock<MockParameterSpace>>()),
      physical_space(std::make_shared<NiceMock<MockPhysicalSpace>>()) {
    b_spline = std::make_unique<spl::BSpline<1>>(physical_space, parameter_space);
  }

 protected:
  std::unique_ptr<spl::BSpline<1>> b_spline;
  std::shared_ptr<NiceMock<MockParameterSpace>> parameter_space;
  std::shared_ptr<NiceMock<MockPhysicalSpace>> physical_space;
};

TEST_F(ABSpline, Returns0_0For0AndDim0) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->Evaluate({ParametricCoordinate{0.0}}, {0})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns0_0For0AndDim1) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->Evaluate({ParametricCoordinate{0.0}}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns4_0For5AndDim0) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->Evaluate({ParametricCoordinate{5.0}}, {0})[0], DoubleEq(4.0));
}

TEST_F(ABSpline, Returns0_0For5AndDim1) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->Evaluate({ParametricCoordinate{5.0}}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns1_5For2_5AndDim0) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->Evaluate({ParametricCoordinate{2.5}}, {0})[0], DoubleEq(1.5));
}

TEST_F(ABSpline, Returns0_0For0_0Dim0AndDer1) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->EvaluateDerivative({ParametricCoordinate{0.0}}, {0}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns1_0For0_0Dim1AndDer1) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->EvaluateDerivative({ParametricCoordinate{0.0}}, {1}, {1})[0], DoubleEq(2.0));
}

TEST_F(ABSpline, Returns12_0For5_0Dim0AndDer1) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->EvaluateDerivative({ParametricCoordinate{5.0}}, {0}, {1})[0], DoubleEq(0.0));
}

TEST_F(ABSpline, Returns0_325For2_25Dim1AndDer1) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->EvaluateDerivative({ParametricCoordinate{2.25}}, {1}, {1})[0], DoubleEq(0.325));
}

TEST_F(ABSpline, CanBeConstructedWithAPhysicalAndAParameterSpace) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THAT(b_spline->Evaluate({ParametricCoordinate{5.0}}, {0})[0], DoubleEq(4.0));
}

TEST_F(ABSpline, ThrowsExceptionForEvaluationAt6_0) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THROW(b_spline->Evaluate({ParametricCoordinate{6.0}}, {0}), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForEvaluationAtMinus1_0) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THROW(b_spline->Evaluate({ParametricCoordinate{-1.0}}, {0}), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForDerivativeEvaluationAt6_0) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THROW(b_spline->EvaluateDerivative({ParametricCoordinate{6.0}}, {0}, {1}), std::runtime_error);
}

TEST_F(ABSpline, ThrowsExceptionForDerivativeEvaluationAtMinus1_0) {  // NOLINT
  mock_parameterSpace(parameter_space);
  mock_physicalSpace(physical_space);
  ASSERT_THROW(b_spline->EvaluateDerivative({ParametricCoordinate{-1.0}}, {0}, {1}), std::runtime_error);
}

class ABSplineWithSplineGenerator : public Test {
 public:
  ABSplineWithSplineGenerator() : degree_{Degree{2}} {
    std::array<baf::KnotVector, 1> knot_vector = {baf::KnotVector(
        {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1},
         ParametricCoordinate{2},
         ParametricCoordinate{3}, ParametricCoordinate{4}, ParametricCoordinate{4}, ParametricCoordinate{5},
         ParametricCoordinate{5},
         ParametricCoordinate{5}})};
    control_points_ = {
        spl::ControlPoint(std::vector<double>({0.0, 0.0})),
        spl::ControlPoint(std::vector<double>({0.0, 1.0})),
        spl::ControlPoint(std::vector<double>({1.0, 1.0})),
        spl::ControlPoint(std::vector<double>({1.5, 1.5})),
        spl::ControlPoint(std::vector<double>({2.0, 1.3})),
        spl::ControlPoint(std::vector<double>({3.0, 2.0})),
        spl::ControlPoint(std::vector<double>({4.0, 1.5})),
        spl::ControlPoint(std::vector<double>({4.0, 0.0}))
    };
    knot_vector_[0] = {std::make_shared<baf::KnotVector>(knot_vector[0])};
    b_spline = std::make_unique<spl::BSpline<1>>(knot_vector_, degree_, control_points_);
  }

 protected:
  std::unique_ptr<spl::BSpline<1>> b_spline;
  baf::KnotVectors<1> knot_vector_;
  std::array<Degree, 1> degree_;
  std::vector<spl::ControlPoint> control_points_;
};

TEST_F(ABSplineWithSplineGenerator, Returns0_0For0AndDim0) {  // NOLINT
  ASSERT_THAT(b_spline->Evaluate({ParametricCoordinate{0.0}}, {0})[0], DoubleEq(0.0));
}
