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

#include <iostream>

#include "gmock/gmock.h"

#include "basis_function_handler.h"
#include "element_integration_point.h"
#include "integration_rule.h"
#include "two_point_gauss_legendre.h"
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

TEST_F(AnIGATestSpline, ElementNURBSBasisFunctions) { // NOLINT
  iga::BasisFunctionHandler basis_function_handler(nurbs_);
  iga::itg::IntegrationRule<2> rule = iga::itg::TwoPointGaussLegendre<2>();
  std::vector<iga::elm::ElementIntegrationPoint> splinelib_element_intg_pnts =
      basis_function_handler.EvaluateAllElementNonZeroNURBSBasisFunctions(0, rule);
  std::vector<double> i1 = {0.240652,0.203999,0.0434425,0.00246914,0.193447,0.163984,0.0349212,0.00198481,0.051834,
                            0.0439395,0.0093571,0.000531828,0.00462963,0.00392451,0.000835743,4.7501e-05};
  std::vector<double> i2 = {0.00462963,0.10015,0.257436,0.128348,0.00372152,0.080505,0.206939,0.103172,0.000997177,
                            0.0215712,0.0554492,0.0276448,8.90643e-05,0.00192667,0.00495252,0.00246914};
  std::vector<double> i3 = {0.00462963,0.00392451,0.000835743,4.7501e-05,0.051834,0.0439395,0.0093571,0.000531828,
                            0.193447,0.163984,0.0349212,0.00198481,0.240652,0.203999,0.0434425,0.00246914};
  std::vector<double> i4 = {8.90643e-05,0.00192667,0.00495252,0.00246914,0.000997177,0.0215712,0.0554492,0.0276448,
                            0.00372152,0.080505,0.206939,0.103172,0.00462963,0.10015,0.257436,0.128348};
  std::vector<std::vector<double>> matlab_element_intg_pnts = {i1, i2, i3, i4};
  for (int i = 0; i < splinelib_element_intg_pnts.size(); ++i) {
    for (int j = 0; j < splinelib_element_intg_pnts[i].GetNonZeroBasisFunctions().size(); ++j) {
      ASSERT_THAT(splinelib_element_intg_pnts[i].GetNonZeroBasisFunctions()[j],
                  DoubleNear(matlab_element_intg_pnts[i][j], 0.00005));
    }
  }
}

TEST_F(AnIGATestSpline, ElementNURBSBasisFunctionDerivatives) { // NOLINT
  iga::BasisFunctionHandler basis_function_handler(nurbs_);
  iga::itg::IntegrationRule<2> rule = iga::itg::TwoPointGaussLegendre<2>();
  std::array<std::vector<iga::elm::ElementIntegrationPoint>, 2> splinelib_element_intg_pnts =
      basis_function_handler.EvaluateAllElementNonZeroNURBSBasisFunctionDerivatives(6, rule);
  std::vector<double> i1_1 = {-0.762835,-0.905235,1.64178,0.0262892,-0.613203,-0.727671,1.31974,0.0211325,-0.164307,
                              -0.194979,0.353624,0.00566243,-0.0146753,-0.0174148,0.0315844,0.000505748};
  std::vector<double> i2_1 = {-0.0547691,-2.03814,1.72675,0.366161,-0.044026,-1.63835,1.38804,0.294338,-0.0117967,
                              -0.438996,0.371925,0.0788675,-0.00105364,-0.0392095,0.033219,0.00704416};
  std::vector<double> i3_1 = {-0.0146753,-0.0174148,0.0315844,0.000505748,-0.164307,-0.194979,0.353624,0.00566243,
                              -0.613203,-0.727671,1.31974,0.0211325,-0.762835,-0.905235,1.64178,0.0262892};
  std::vector<double> i4_1 = {-0.00105364,-0.0392095,0.033219,0.00704416,-0.0117967,-0.438996,0.371925,0.0788675,
                              -0.044026,-1.63835,1.38804,0.294338,-0.0547691,-2.03814,1.72675,0.366161};
  std::vector<double> i1_2 = {-0.152567,-3.03375,-0.544322,-0.00140883,0.0708066,1.40797,0.252621,0.000653841,
                              0.0708066,1.40797,0.252621,0.000653841,0.0109538,0.217814,0.0390805,0.00010115};
  std::vector<double> i2_2 = {-0.00293507,-2.31552,-1.34036,-0.0732322,0.00136217,1.07464,0.622065,0.0339872,
                              0.00136217,1.07464,0.622065,0.0339872,0.000210728,0.166247,0.0962338,0.00525783};
  std::vector<double> i3_2 = {-0.0109538,-0.217814,-0.0390805,-0.00010115,-0.0708066,-1.40797,-0.252621,-0.000653841,
                              -0.0708066,-1.40797,-0.252621,-0.000653841,0.152567,3.03375,0.544322,0.00140883};
  std::vector<double> i4_2 = {-0.000210728,-0.166247,-0.0962338,-0.00525783,-0.00136217,-1.07464,-0.622065,-0.0339872,
                              -0.00136217,-1.07464,-0.622065,-0.0339872,0.00293507,2.31552,1.34036,0.0732322};
  std::vector<std::vector<std::vector<double>>> matlab_element_intg_pnts = {{i1_1, i2_1, i3_1, i4_1},
                                                                            {i1_2, i2_2, i3_2, i4_2}};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < splinelib_element_intg_pnts[i].size(); ++j) {
      for (int k = 0; k < splinelib_element_intg_pnts[i][j].GetNonZeroBasisFunctions().size(); ++k) {
        ASSERT_THAT(splinelib_element_intg_pnts[i][j].GetNonZeroBasisFunctions()[k],
                  DoubleNear(matlab_element_intg_pnts[i][j][k], 0.00005));
      }
    }
  }
}
