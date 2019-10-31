/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "gmock/gmock.h"

#include "src/spl/b_spline.h"
#include "src/spl/nurbs.h"
#include "nurbs_3d_mocking.h"
#include "src/util/numeric_settings.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

using namespace splinelib::src;

/* 3-dimensional nurbs spline with following properties :
 * KnotVector = {{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}, {0, 0, 1, 1}}
 * ControlPoints = {{0, 0}, {1, 0}, {3, 0}, {-1, 0.5}, {2, 2}, {4, 1},
 *                  {0, 2}, {-1, 0.5}, {2, 2}, {4, 1}, {0, 2}, {5, 2}}
 * Weights = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
*/

class A3DNurbsWithAllWeights1 : public Test {
 public:
  A3DNurbsWithAllWeights1() :
      parameter_space(std::make_shared<NiceMock<MockParameterSpace3d>>()),
      w_physical_space(std::make_shared<NiceMock<MockWeightedPhysicalSpace3d>>()),
      physical_space(std::make_shared<NiceMock<MockPhysicalSpace3d>>()) {
    spl::NURBSGenerator<3> nurbs_generator(w_physical_space, parameter_space);
    spl::BSplineGenerator<3> bspline_generator(physical_space, parameter_space);

    nurbs_ = std::make_unique<spl::NURBS<3>>(nurbs_generator);
    bspline_ = std::make_unique<spl::BSpline<3>>(bspline_generator);
  }

 protected:
  std::unique_ptr<spl::NURBS<3>> nurbs_;
  std::unique_ptr<spl::BSpline<3>> bspline_;
  std::shared_ptr<NiceMock<MockParameterSpace3d>> parameter_space;
  std::shared_ptr<NiceMock<MockWeightedPhysicalSpace3d>> w_physical_space;
  std::shared_ptr<NiceMock<MockPhysicalSpace3d>> physical_space;
};

TEST_F(A3DNurbsWithAllWeights1, ReturnsSameDerivativeAs3DBSplineFor0_5And0_5And0_5AndDerivatives1And1And0) { // NOLINT
  mock_parameterSpace_nurbs3d(parameter_space);
  mock_physicalSpace3d(physical_space);
  mock_weightedPhysicalSpace3d(w_physical_space);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParametricCoordinate{0.5}, ParametricCoordinate{0.5},
                                          ParametricCoordinate{0.5}}, {0},
                                         {1, 1, 0})[0],
              DoubleEq(bspline_->EvaluateDerivative({ParametricCoordinate{0.5}, ParametricCoordinate{0.5},
                                                     ParametricCoordinate{0.5}}, {0}, {1, 1, 0})[0]));
}

TEST_F(A3DNurbsWithAllWeights1, ReturnsSameDerivativeAs3DBSplineFor0_5And0_8And0_1AndDerivatives1And1And1) { // NOLINT
  mock_parameterSpace_nurbs3d(parameter_space);
  mock_physicalSpace3d(physical_space);
  mock_weightedPhysicalSpace3d(w_physical_space);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParametricCoordinate{0.5}, ParametricCoordinate{0.8},
                                          ParametricCoordinate{0.1}}, {0},
                                         {1, 1, 1})[0],
              DoubleEq(bspline_->EvaluateDerivative({ParametricCoordinate{0.5}, ParametricCoordinate{0.8},
                                                     ParametricCoordinate{0.1}}, {0}, {1, 1, 1})[0]));
}

TEST_F(A3DNurbsWithAllWeights1, ReturnsSameDerivativeAs3DBSplineFor0_5And0_8And0_1AndDerivatives1And2And1) { // NOLINT
  mock_parameterSpace_nurbs3d(parameter_space);
  mock_physicalSpace3d(physical_space);
  mock_weightedPhysicalSpace3d(w_physical_space);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParametricCoordinate{0.5}, ParametricCoordinate{0.8},
                                          ParametricCoordinate{0.1}}, {0},
                                         {1, 2, 1})[0],
              DoubleNear(bspline_->EvaluateDerivative({ParametricCoordinate{0.5}, ParametricCoordinate{0.8},
                                                       ParametricCoordinate{0.1}}, {0}, {1, 2, 1})[0],
                         util::numeric_settings::GetEpsilon<double>()));
}
