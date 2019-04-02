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

#include "basis_function_handler.h"
#include "two_point_gauss_legendre.h"
#include "nurbs.h"

const static std::array<baf::KnotVector, 2> knot_vector =
      {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.4}, ParamCoord{0.5},
                        ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}),
       baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{0.5},
                        ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})};

const static std::array<Degree, 2> degree = {Degree{3}, Degree{3}};

const static std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1};

const static std::vector<baf::ControlPoint> control_points = {
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
    baf::ControlPoint(std::vector<double>({-1.0, 1.0, 0.0}))};

const static std::array<std::shared_ptr<baf::KnotVector>, 2> kv_ptr = {
    std::make_shared<baf::KnotVector>(knot_vector[0]), std::make_shared<baf::KnotVector>(knot_vector[1])};
const static std::shared_ptr<spl::NURBS<2>> nurbs_ = std::make_shared<spl::NURBS<2>>(kv_ptr, degree, control_points,
    weights);
const static iga::itg::IntegrationRule rule = iga::itg::TwoPointGaussLegendre();
const static iga::BasisFunctionHandler<2> basis_function_handler = iga::BasisFunctionHandler<2>(nurbs_);

int main() {
  std::vector<iga::elm::ElementIntegrationPoint<2>> splinelib_element_intg_pnts =
      basis_function_handler.EvaluateAllElementNonZeroNURBSBasisFunctions(0, rule);
  std::vector<double> i1 = {0.240652, 0.203999, 0.0434425, 0.00246914, 0.193447, 0.163984, 0.0349212, 0.00198481,
                            0.051834, 0.0439395, 0.0093571, 0.0005318, 0.00462963, 0.00392451, 0.000835743, 4.7501e-05};
  std::vector<double> i2 = {0.004629, 0.10015, 0.257436, 0.128348, 0.0037215, 0.080505, 0.206939, 0.103172, 0.000997177,
                            0.0215712, 0.0554492, 0.0276448, 8.90643e-05, 0.00192667, 0.00495252, 0.00246914};
  std::vector<double> i3 = {0.0046296, 0.0039245, 0.000835743, 4.7501e-05, 0.051834, 0.0439395, 0.0093571, 0.000531828,
                            0.193447, 0.163984, 0.0349212, 0.00198481, 0.240652, 0.203999, 0.0434425, 0.00246914};
  std::vector<double> i4 = {8.90643e-05, 0.00192667, 0.00495252, 0.00246914, 0.0009972, 0.0215712, 0.0554492, 0.0276448,
                            0.00372152, 0.080505, 0.206939, 0.103172, 0.00462963, 0.10015, 0.257436, 0.128348};
  std::vector<std::vector<double>> matlab_elm_intg_pnts = {i1, i2, i3, i4};
  for (u_int64_t i = 0; i < splinelib_element_intg_pnts.size(); ++i) {
    for (uint64_t j = 0; j < splinelib_element_intg_pnts[i].GetNonZeroBasisFunctions().size(); ++j) {
      if (abs(splinelib_element_intg_pnts[i].GetNonZeroBasisFunctions()[j] -  matlab_elm_intg_pnts[i][j]) > 0.00005) {
        std::cout << std::endl << "ERROR!" << std::endl;
      }
    }
  }

  return 0;
}
