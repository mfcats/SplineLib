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
#include "b_spline_2d_mocking.h"

using testing::Test;
using ::testing::NiceMock;

using namespace splinelib::src;

/* 2-dimensional spline with following properties :
 * KnotVector1 = {0, 0, 0, 1, 1, 1}
 * KnotVector2 = {0, 0, 0, 1, 1, 1}
 * ControlPoints = {{-1, -1, 0}, {0, -1, 0}, {1, -1, 0}, {-1, 0, 0}, {0, 0, 1},
 *                  {1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {1, 1, 0}}
*/

class A2DBSpline : public Test {
 public:
  A2DBSpline() :
      parameter_space(std::make_shared<NiceMock<Mock2dParameterSpace>>()),
      physical_space(std::make_shared<NiceMock<Mock2dPhysicalSpace>>()) {
    spl::BSplineGenerator<2> b_spline_generator(physical_space, parameter_space);
    b_spline = std::make_unique<spl::BSpline<2>>(b_spline_generator);
  }

 protected:
  std::unique_ptr<spl::BSpline<2>> b_spline;
  std::shared_ptr<NiceMock<Mock2dParameterSpace>> parameter_space;
  std::shared_ptr<NiceMock<Mock2dPhysicalSpace>> physical_space;
};

TEST_F(A2DBSpline, Corner) {  // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {0})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {1})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {2})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, EdgeDim0) {  // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.33333}}, {0})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.33333}}, {1})[0],
              -0.33333,
              0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.33333}}, {2})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, EdgeDim1) {  // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.33333}, ParametricCoordinate{0.0}}, {0})[0],
              -0.33333,
              0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.33333}, ParametricCoordinate{0.0}}, {1})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.33333}, ParametricCoordinate{0.0}}, {2})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, Center) {  // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.5}, ParametricCoordinate{0.5}}, {0})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.5}, ParametricCoordinate{0.5}}, {1})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.5}, ParametricCoordinate{0.5}}, {2})[0], 0.25, 0.00005);
}

TEST_F(A2DBSpline, Random) {  // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.75}, ParametricCoordinate{0.25}}, {0})[0], 0.5, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.75}, ParametricCoordinate{0.25}}, {1})[0], -0.5, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParametricCoordinate{0.75}, ParametricCoordinate{0.25}}, {2})[0], 0.14063, 0.00005);
}

TEST_F(A2DBSpline, CornerDer10) {  // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {0}, {1, 0})[0],
              2.0,
              0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {1}, {1, 0})[0],
              0.0,
              0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {2}, {1, 0})[0],
              0.0,
              0.00005);
}

TEST_F(A2DBSpline, CornerDer01) {  // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {0}, {0, 1})[0],
              0.0,
              0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {1}, {0, 1})[0],
              2.0,
              0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {2}, {0, 1})[0],
              0.0,
              0.00005);
}

TEST_F(A2DBSpline, EdgeDim0Der10) {  // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.0}, ParametricCoordinate{0.33333}}, {0}, {1, 0})[0], 2.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.0}, ParametricCoordinate{0.33333}}, {1}, {1, 0})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.0}, ParametricCoordinate{0.33333}}, {2}, {1, 0})[0], 0.888889, 0.00005);
}

TEST_F(A2DBSpline, EdgeDim0Der01) {  // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.0}, ParametricCoordinate{0.33333}}, {0}, {0, 1})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.0}, ParametricCoordinate{0.33333}}, {1}, {0, 1})[0], 2.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.0}, ParametricCoordinate{0.33333}}, {2}, {0, 1})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, CenterDer10) {  // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParametricCoordinate{0.5}, ParametricCoordinate{0.5}}, {0}, {1, 0})[0],
              2.0,
              0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParametricCoordinate{0.5}, ParametricCoordinate{0.5}}, {1}, {1, 0})[0],
              0.0,
              0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParametricCoordinate{0.5}, ParametricCoordinate{0.5}}, {2}, {1, 0})[0],
              0.0,
              0.00005);
}

TEST_F(A2DBSpline, RandomDer10) {  // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.75}, ParametricCoordinate{0.25}}, {0}, {1, 0})[0], 2.000, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.75}, ParametricCoordinate{0.25}}, {1}, {1, 0})[0], 0.000, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.75}, ParametricCoordinate{0.25}}, {2}, {1, 0})[0], -0.375, 0.00005);
}

TEST_F(A2DBSpline, RandomDer01) {  // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.75}, ParametricCoordinate{0.25}}, {0}, {0, 1})[0], 0.000, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.75}, ParametricCoordinate{0.25}}, {1}, {0, 1})[0], 2.000, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.75}, ParametricCoordinate{0.25}}, {2}, {0, 1})[0], 0.375, 0.00005);
}

TEST_F(A2DBSpline, RandomDer12) {  // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.75}, ParametricCoordinate{0.25}}, {0}, {1, 2})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.75}, ParametricCoordinate{0.25}}, {1}, {1, 2})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative(
      {ParametricCoordinate{0.75}, ParametricCoordinate{0.25}}, {2}, {1, 2})[0], 4.0, 0.00005);
}
