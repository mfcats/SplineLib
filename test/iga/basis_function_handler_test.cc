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
#include <vector>

#include "gmock/gmock.h"

#include "basis_function_handler.h"
#include "nurbs.h"
#include "test_spline.h"

using testing::DoubleNear;

TEST_F(AnIGATestSpline, TestBasisFunctionHandler) { // NOLINT
  iga::BasisFunctionHandler basis_function_handler(nurbs_);
  std::vector<double> splinelib_nurbs_baf =
      basis_function_handler.EvaluateAllNonZeroNURBSBasisFunctions(std::array<ParamCoord, 2>({ParamCoord{0.55},
                                                                                              ParamCoord{0.3}}));
  std::vector<double> matlab_nurbs_baf = {0.000667, 0.046933, 0.01608, 0.00032, 0.003, 0.2112, 0.07236, 0.00144, 0.0045,
                                          0.3168, 0.10854, 0.00216, 0.00225, 0.1584, 0.05427, 0.00108};
  for (int i = 0; i < splinelib_nurbs_baf.size(); ++i) {
    ASSERT_THAT(splinelib_nurbs_baf[i], DoubleNear(matlab_nurbs_baf[i], 0.00005));
  }
}

TEST_F(AnIGATestSpline, TestBasisFunctionHandlerDerivativeParam) { // NOLINT
  iga::BasisFunctionHandler basis_function_handler(nurbs_);
  std::array<std::vector<double>, 2> splinelib_nurbs_baf_der_param =
      basis_function_handler.EvaluateAllNonZeroNURBSBasisFunctionDerivatives(
          std::array<ParamCoord, 2>({ParamCoord{0.7364}, ParamCoord{0.3892}}));
  std::array<std::vector<double>, 2> matlab_nurbs_baf_der_param =
      std::array<std::vector<double>, 2>({
        std::vector<double>({-0.018903, -0.016112, 0.025525, 0.009490, -0.199202, -0.169791, 0.268985, 0.100008,
                             -0.699725, -0.596412, 0.944847, 0.351291, -0.819293, -0.698326, 1.106301, 0.411319}),
        std::vector<double>({-0.044972, -0.140697, -0.097287, -0.011683, -0.270971, -0.847739, -0.586179, -0.070393,
                             -0.238953, -0.747572, -0.516917, -0.062076, 0.554896, 1.736009, 1.200382, 0.144152})});
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < splinelib_nurbs_baf_der_param.at(i).size(); ++j) {
    ASSERT_THAT(splinelib_nurbs_baf_der_param.at(i).at(j), DoubleNear(matlab_nurbs_baf_der_param.at(i).at(j), 0.00005));
    }
  }
}

TEST_F(AnIGATestSpline, TestBasisFunctionHandlerDerivativePhysical) { // NOLINT
  iga::BasisFunctionHandler basis_function_handler(nurbs_);
  std::array<std::vector<double>, 2> splinelib_nurbs_baf_der_phy =
      basis_function_handler.GetDrDx(std::array<ParamCoord, 2>({ParamCoord{0.367}, ParamCoord{0.893}}));
  std::array<std::vector<double>, 2> matlab_nurbs_baf_der_phy =
      std::array<std::vector<double>, 2>({
        std::vector<double>({7.281165e-05, 0.011839, 0.064344, 0.053414, 0.000462, 0.075125, 0.408316, 0.338955,
                             0.000447, 0.072741, 0.395356, 0.328197, -0.000982, -0.159705, -0.868017, -0.720566}),
        std::vector<double>({-0.000587, -0.021549, -0.017258, 0.038123, -0.006460, -0.236904, -0.187217, 0.422517,
                             -0.023715, -0.868134, -0.676805, 1.560847, -0.02902, -1.060418, -0.815362, 1.921943})});
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < splinelib_nurbs_baf_der_phy.at(i).size(); ++j) {
      ASSERT_THAT(splinelib_nurbs_baf_der_phy.at(i).at(j), DoubleNear(matlab_nurbs_baf_der_phy.at(i).at(j), 0.00005));
    }
  }
}
