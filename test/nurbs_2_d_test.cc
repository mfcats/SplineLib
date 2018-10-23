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

#include "nurbs.h"
#include "b_spline.h"
#include "2d_nurbs_mocking.h"

using testing::Test;
using ::testing::NiceMock;
using testing::DoubleEq;
using testing::DoubleNear;

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


/* 2-dimensional nurbs spline with following properties :
 * KnotVector = {{0, 0, 0, 2, 2, 2}, {0, 0, 0, 2, 2, 2}}
 * ControlPoints = {{1, 2}, {2, 2}, {4, 2}, {0, 2.5}, {3, 4}, {5, 3}, {1, 4}, {3.5, 5.5}, {6, 4}}
 * Weights = {1, 1, 1, 1, 1, 1, 1, 1, 1}
*/

class A2DNurbsWithAllWeights1 : public Test {
 public:
  A2DNurbsWithAllWeights1() :
      parameter_space_m(std::make_shared<NiceMock<MockParameterSpace2>>()),
      w_physical_space_m(std::make_shared<NiceMock<MockWeightedPhysicalSpace2>>()),
      physical_space_m(std::make_shared<NiceMock<MockPhysicalSpace2>>()) {
    spl::NURBSGenerator<2> nurbs_generator(w_physical_space_m, parameter_space_m);
    nurbs_ = std::make_unique<spl::NURBS<2>>(nurbs_generator);
    spl::BSplineGenerator<2> bspline_generator(physical_space_m, parameter_space_m);
    bspline_ = std::make_unique<spl::BSpline<2>>(bspline_generator);
  }

 protected:
  std::unique_ptr<spl::NURBS<2>> nurbs_;
  std::unique_ptr<spl::BSpline<2>> bspline_;
  std::shared_ptr<NiceMock<MockParameterSpace2>> parameter_space_m;
  std::shared_ptr<NiceMock<MockWeightedPhysicalSpace2>>w_physical_space_m;
  std::shared_ptr<NiceMock<MockPhysicalSpace2>> physical_space_m;
};

TEST_F(A2DNurbsWithAllWeights1, ReturnsSameDerivativeAs2DBSplineFor0_5And0_5AndDerivatives1And1) { // NOLINT
  mock_weightedPhysicalSpace(w_physical_space_m);
  mock_parameterSpace_nurbs(parameter_space_m);
  mock_physicalSpace(physical_space_m);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}}, {0}, {1, 1})[0],
              DoubleNear(bspline_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}}, {0}, {1, 1})[0], 0.000001));

}

TEST_F(A2DNurbsWithAllWeights1, ReturnsSameDerivativeAs2DBSplineFor0_0And0_7AndDerivatives1And1) { // NOLINT
  mock_weightedPhysicalSpace(w_physical_space_m);
  mock_parameterSpace_nurbs(parameter_space_m);
  mock_physicalSpace(physical_space_m);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.7}}, {0}, {1, 1})[0],
              DoubleNear(bspline_->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.7}}, {0}, {1, 1})[0], 0.000001));
}

TEST_F(A2DNurbsWithAllWeights1, ReturnsSameDerivativeAs2DBSplineFor0_0And0_7AndDerivatives2And1) { // NOLINT
  mock_weightedPhysicalSpace(w_physical_space_m);
  mock_parameterSpace_nurbs(parameter_space_m);
  mock_physicalSpace(physical_space_m);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.7}}, {0}, {2, 1})[0],
              DoubleNear(bspline_->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.7}}, {0}, {2, 1})[0], 0.000001));
}
