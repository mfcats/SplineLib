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

#include "gmock/gmock.h"
#include "nurbs.h"
#include "basis_function_handler.h"

using testing::DoubleNear;
using testing::Test;

class ABasisFunctionHandler : public Test {
 public:
  std::array<baf::KnotVector, 2> knot_vector =
      {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.4}, ParamCoord{0.5},
                        ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}),
       baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{0.5},
                        ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})};

  std::array<Degree, 2> degree = {Degree{3}, Degree{3}};

  std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1};

  std::vector<baf::ControlPoint> control_points = {
      baf::ControlPoint(std::vector<double>({0.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.16, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.32, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.48, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.64, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.8, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -1.0, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.16, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.32, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.48, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.64, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.8, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -0.66, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.16, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.32, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.48, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.64, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.8, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -0.33, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.16, 0.16, 0.0})),
      baf::ControlPoint(std::vector<double>({0.32, 0.32, 0.0})),
      baf::ControlPoint(std::vector<double>({0.48, 0.48, 0.0})),
      baf::ControlPoint(std::vector<double>({0.64, 0.64, 0.0})),
      baf::ControlPoint(std::vector<double>({0.8, 0.8, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0})),

      baf::ControlPoint(std::vector<double>({-0.33, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.16, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.32, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.48, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.64, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.8, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 1.0, 0.0})),

      baf::ControlPoint(std::vector<double>({-0.66, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.16, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.32, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.48, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.64, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.8, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 1.0, 0.0})),

      baf::ControlPoint(std::vector<double>({-1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.16, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.32, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.48, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.64, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.8, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 1.0, 0.0}))
  };

  std::array<std::shared_ptr<baf::KnotVector>, 2> kv_ptr = {std::make_shared<baf::KnotVector>(knot_vector[0]),
                                                            std::make_shared<baf::KnotVector>(knot_vector[1])};
  std::shared_ptr<spl::NURBS<2>> nurbs_ = std::make_shared<spl::NURBS<2>>(kv_ptr, degree, control_points, weights);
};

TEST_F(ABasisFunctionHandler, TestBasisFunctionHandler) { // NOLINT
  iga::BasisFunctionHandler basis_function_handler(nurbs_);
  std::vector<double> splinelib_nurbs_baf =
      basis_function_handler.EvaluateAllNonZeroNURBSBasisFunctions(std::array<ParamCoord,2>({ParamCoord{0.55},
                                                                                             ParamCoord{0.3}}));
  std::vector<double> matlab_nurbs_baf = {0.000667, 0.046933, 0.01608, 0.00032, 0.003, 0.2112, 0.07236, 0.00144, 0.0045,
                                          0.3168, 0.10854, 0.00216, 0.00225, 0.1584, 0.05427, 0.00108};
  for (int i = 0; i < splinelib_nurbs_baf.size(); ++i) {
    ASSERT_THAT(splinelib_nurbs_baf[i], DoubleNear(matlab_nurbs_baf[i], 0.00005));
  }
}

