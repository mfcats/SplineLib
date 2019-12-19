/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/spl/nurbs.h"
#include "src/util/numeric_settings.h"
#include "test/spl/mocking/nurbs_1d_mocking.h"

using testing::Test;
using testing::DoubleNear;

using namespace splinelib::src;

/* 1-dimensional nurbs spline with following properties :
 * KnotVector = {0, 0, 0, 1, 2, 3, 3, 3}
 * ControlPoints = {{0, 0}, {1, 1}, {3, 2}, {4, 1}, {5, -1}}
 * Weights = {1, 4, 1, 1, 1}
*/
class NurbsEx4_1 : public Test {
 public:
  NurbsEx4_1() :
      parameter_space(std::make_shared<NiceMock<MockParameterSpace14111>>()),
      w_physical_space(std::make_shared<NiceMock<MockWeightedPhysicalSpace14111>>()) {
    nurbs = std::make_unique<spl::NURBS<1>>(w_physical_space, parameter_space);
  }

 protected:
  std::unique_ptr<spl::NURBS<1>> nurbs;
  std::shared_ptr<NiceMock<MockParameterSpace14111>> parameter_space;
  std::shared_ptr<NiceMock<MockWeightedPhysicalSpace14111>> w_physical_space;
};

TEST_F(NurbsEx4_1, Returns1_4For1AndDim0) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{1.0}}, {0})[0],
              DoubleNear(1.4, util::numeric_settings::GetEpsilon<double>()));
}

TEST_F(NurbsEx4_1, Returns1_2For1AndDim1) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{1.0}}, {1})[0],
              DoubleNear(1.2, util::numeric_settings::GetEpsilon<double>()));
}

TEST_F(NurbsEx4_1, Returns1For1AndWeight) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{1.0}}, {2})[0],
              DoubleNear(2.5, util::numeric_settings::GetEpsilon<double>()));
}

/* 1-dimensional nurbs spline with following properties :
 * KnotVector = {0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 0.95, 1.0, 1.0, 1.0}
 * ControlPoints = {{0,5, 3.0, 1.0}, {1.5, 5.5, 4.0}, {4.5, 5.5, 0.1}, {3.0, 1.5, 2.0}
 *                  {7.5, 1.5, 3.5}, {6.0, 4.0, 5.3}, {8.5, 4.5, 0.0}}
 * Weights = {1.0, 0.9, 0.7, 0.5, 0.8, 1.2, 2.0}
*/
class ANurbs : public Test {
 public:
  ANurbs() :
      parameter_space(std::make_shared<NiceMock<MockParameterSpace1009>>()),
      w_physical_space(std::make_shared<NiceMock<MockWeightedPhysicalSpace1009>>()) {
    nurbs = std::make_unique<spl::NURBS<1>>(w_physical_space, parameter_space);
  }

 protected:
  std::unique_ptr<spl::NURBS<1>> nurbs;
  std::shared_ptr<NiceMock<MockParameterSpace1009>> parameter_space;
  std::shared_ptr<NiceMock<MockWeightedPhysicalSpace1009>> w_physical_space;
};

TEST_F(ANurbs, ReturnsCorrectCurvePointForFirstKnot) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{0.0}}, {0})[0],
              DoubleNear(0.5, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{0.0}}, {1})[0],
              DoubleNear(3.0, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{0.0}}, {2})[0],
              DoubleNear(1.0, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{0.0}}, {3})[0],
              DoubleNear(1.0, util::numeric_settings::GetEpsilon<double>()));
}

TEST_F(ANurbs, ReturnsCorrectCurvePointForInnerKnot) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{0.25}}, {0})[0],
              DoubleNear(2.8125, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{0.25}}, {1})[0],
              DoubleNear(5.5, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{0.25}}, {2})[0],
              DoubleNear(2.29375, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{0.25}}, {3})[0],
              DoubleNear(0.8, util::numeric_settings::GetEpsilon<double>()));
}

TEST_F(ANurbs, ReturnsCorrectCurvePointForValueBetweenTwoKnots) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{1.0 / 3.0}}, {0})[0],
              DoubleNear(3.625, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{1.0 / 3.0}}, {1})[0], DoubleNear(5.34848, 0.000005));
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{1.0 / 3.0}}, {2})[0], DoubleNear(1.23561, 0.000005));
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{1.0 / 3.0}}, {3})[0], DoubleNear(0.73333, 0.000005));
}

TEST_F(ANurbs, ReturnsCorrectCurvePointForLastKnot) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{1.0}}, {0})[0],
              DoubleNear(8.5, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{1.0}}, {1})[0],
              DoubleNear(4.5, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{1.0}}, {2})[0],
              DoubleNear(0.0, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{1.0}}, {3})[0],
              DoubleNear(2.0, util::numeric_settings::GetEpsilon<double>()));
}

TEST_F(ANurbs, ThrowsExceptionForEvaluationAt1_2) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  ASSERT_THROW(nurbs->Evaluate({ParametricCoordinate{1.2}}, {0}), std::runtime_error);
}

TEST_F(ANurbs, ThrowsExceptionForEvaluationAtMinus0_1) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  ASSERT_THROW(nurbs->Evaluate({ParametricCoordinate{-0.1}}, {0}), std::runtime_error);
}

/* 1-dimensional nurbs spline with following properties :
 * KnotVector = {0, 0, 0, 1, 1, 1}
 * ControlPoints = {{1, 0}, {1. 1}, {0, 1}}
 * Weights = {1, 1, 2}
*/
class NurbsDerivativeEx4_2 : public Test {
 public:
  NurbsDerivativeEx4_2() :
      parameter_space(std::make_shared<NiceMock<MockParameterSpace112>>()),
      w_physical_space(std::make_shared<NiceMock<MockWeightedPhysicalSpace112>>()) {
    nurbs = std::make_unique<spl::NURBS<1>>(w_physical_space, parameter_space);
  }

 protected:
  std::unique_ptr<spl::NURBS<1>> nurbs;
  std::shared_ptr<NiceMock<MockParameterSpace112>> parameter_space;
  std::shared_ptr<NiceMock<MockWeightedPhysicalSpace112>> w_physical_space;
};

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForFirstDerivativeAtFirstKnot) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{0.0}}, {0}, {1})[0], 0.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{0.0}}, {1}, {1})[0], 2.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForFirstDerivativeAtValueBetweenKnots) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{0.5}}, {0}, {1})[0], -1.28);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{0.5}}, {1}, {1})[0], 0.96);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForFirstDerivativeAtLastKnot) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{1.0}}, {0}, {1})[0], -1.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{1.0}}, {1}, {1})[0], 0.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForSecondDerivativeAtFirstKnot) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{0.0}}, {0}, {2})[0], -4.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{0.0}}, {1}, {2})[0], 0.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForSecondDerivativeAtValueBetweenKnots) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{0.5}}, {0}, {2})[0],
              DoubleNear(-0.512, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{0.5}}, {1}, {2})[0],
              DoubleNear(-2.816, util::numeric_settings::GetEpsilon<double>()));
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForSecondDerivativeAtLastKnot) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{1.0}}, {0}, {2})[0], 1.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{1.0}}, {1}, {2})[0], -1.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForThirdDerivativeAtFirstKnot) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{0.0}}, {0}, {3})[0], 0.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{0.0}}, {1}, {3})[0], -12.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForThirdDerivativeAtValueBetweenKnots) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{0.5}}, {0}, {3})[0],
              DoubleNear(7.3728, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{0.5}}, {1}, {3})[0],
              DoubleNear(2.1504, util::numeric_settings::GetEpsilon<double>()));
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForThirdDerivativeAtLastKnot) {  // NOLINT
  mock_parameterSpace_nurbs(parameter_space);
  mock_weightedPhysicalSpace(w_physical_space);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{1.0}}, {0}, {3})[0], 0.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParametricCoordinate{1.0}}, {1}, {3})[0], 3.0);
}

class ANURBSWithSplineGenerator : public Test {
 public:
  ANURBSWithSplineGenerator() {
    baf::KnotVectors<1> knot_vector =
        {std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1},
             ParametricCoordinate{2},
             ParametricCoordinate{3}, ParametricCoordinate{3}, ParametricCoordinate{3}}))};
    std::array<Degree, 1> degree = {Degree{2}};
    std::vector<double> weights = {1, 4, 1, 1, 1};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint(std::vector<double>({0.0, 0.0})),
        spl::ControlPoint(std::vector<double>({1.0, 1.0})),
        spl::ControlPoint(std::vector<double>({3.0, 2.0})),
        spl::ControlPoint(std::vector<double>({4.0, 1.0})),
        spl::ControlPoint(std::vector<double>({5.0, -1.0}))
    };
    nurbs = std::make_unique<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::unique_ptr<spl::NURBS<1>> nurbs;
};

TEST_F(ANURBSWithSplineGenerator, Returns1_4For1AndDim0) {  // NOLINT
  ASSERT_THAT(nurbs->Evaluate({ParametricCoordinate{1.0}}, {0})[0],
              DoubleNear(1.4, util::numeric_settings::GetEpsilon<double>()));
}
