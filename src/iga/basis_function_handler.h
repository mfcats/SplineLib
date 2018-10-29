/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_BASIS_FUNCTION_HANDLER_H_
#define SRC_IGA_BASIS_FUNCTION_HANDLER_H_

#include <math.h>
#include <array>
#include <vector>

#include "element_integration_point.h"
#include "mapping_handler.h"
#include "nurbs.h"

namespace iga {
class BasisFunctionHandler {
 public:
  explicit BasisFunctionHandler(std::shared_ptr<spl::NURBS<2>> spl) : spline_(std::move(spl)) {}

  std::vector<double> EvaluateAllNonZeroNURBSBasisFunctions(std::array<ParamCoord, 2> param_coord) {
    std::array<std::vector<double>, 2> basis_functions;
    basis_functions[0] = spline_->EvaluateAllNonZeroBasisFunctions(0, param_coord[0]);
    basis_functions[1] = spline_->EvaluateAllNonZeroBasisFunctions(1, param_coord[1]);
    std::vector<double> nurbs_basis_functions;
    double sum = 0;
    for (int i = 0; i < basis_functions[1].size(); ++i) {
      for (int j = 0; j < basis_functions[0].size(); ++j) {
        double temp = basis_functions[0][j] * basis_functions[1][i] * spline_->GetWeight({j, i});
        sum += temp;
        nurbs_basis_functions.emplace_back(temp);
      }
    }
    for (int i = 0; i < nurbs_basis_functions.size(); ++i) {
      nurbs_basis_functions[i] = nurbs_basis_functions[i] / sum;
    }
    return nurbs_basis_functions;
  }

  std::array<std::vector<double>, 2> EvaluateAllNonZeroNURBSBasisFunctionDerivatives(
      std::array<ParamCoord, 2> param_coord) {
    std::vector<double> nurbs_basis_functions;
    std::array<std::vector<double>, 2> nurbs_basis_function_derivatives;
    std::array<std::vector<double>, 2> basis_functions;
    std::array<std::vector<double>, 2> basis_function_derivatives;
    basis_functions[0] = spline_->EvaluateAllNonZeroBasisFunctions(0, param_coord[0]);
    basis_functions[1] = spline_->EvaluateAllNonZeroBasisFunctions(1, param_coord[1]);
    basis_function_derivatives[0] = spline_->EvaluateAllNonZeroBasisFunctionDerivatives(0, param_coord[0], 1);
    basis_function_derivatives[1] = spline_->EvaluateAllNonZeroBasisFunctionDerivatives(1, param_coord[1], 1);
    double sum_baf = 0;
    double sum_der_xi = 0;
    double sum_der_eta = 0;
    for (int i = 0; i < basis_function_derivatives[1].size(); ++i) {
      for (int j = 0; j < basis_function_derivatives[0].size(); ++j) {
        double temp = basis_functions[0][j] * basis_functions[1][i] * spline_->GetWeight({j, i});
        double temp_der_xi = basis_function_derivatives[0][j] * basis_functions[1][i] * spline_->GetWeight({j, i});
        double temp_der_eta = basis_functions[0][j] * basis_function_derivatives[1][i] * spline_->GetWeight({j, i});
        sum_baf += temp;
        sum_der_xi += temp_der_xi;
        sum_der_eta += temp_der_eta;
        nurbs_basis_functions.emplace_back(temp);
        nurbs_basis_function_derivatives[0].emplace_back(temp_der_xi);
        nurbs_basis_function_derivatives[1].emplace_back(temp_der_eta);
      }
    }
    for (int i = 0; i < basis_functions.size(); ++i) {
      nurbs_basis_functions[i] = nurbs_basis_functions[i] / sum_baf;
      nurbs_basis_function_derivatives[0][i] = (nurbs_basis_function_derivatives[0][i] * sum_baf -
          nurbs_basis_functions[i] * sum_der_xi) / pow(sum_baf, 2);
      nurbs_basis_function_derivatives[1][i] = (nurbs_basis_function_derivatives[1][i] * sum_baf -
          nurbs_basis_functions[i] * sum_der_eta) / pow(sum_baf, 2);
    }
    return nurbs_basis_function_derivatives;
  }

  std::array<std::vector<double>, 2> GetDrDx(std::array<ParamCoord, 2> param_coord) {
    std::array<std::vector<double>, 2> dr_dx;
    MappingHandler mapping_handler(spline_);
    std::array<std::vector<double>, 2> dr_dxi = EvaluateAllNonZeroNURBSBasisFunctionDerivatives(param_coord);
    for (int i = 0; i < dr_dxi[0].size(); ++i) {
      dr_dx[0].emplace_back(dr_dxi[0][i] * mapping_handler.GetDxiDx(param_coord)[0][0]
                            + dr_dxi[1][i] * mapping_handler.GetDxiDx(param_coord)[1][0]);
      dr_dx[1].emplace_back(dr_dxi[0][i] * mapping_handler.GetDxiDx(param_coord)[0][1]
                            + dr_dxi[1][i] * mapping_handler.GetDxiDx(param_coord)[1][1]);
    }
    return dr_dx;
  }

 private:
  std::shared_ptr<spl::NURBS<2>> spline_;
};
}  // namespace iga

#endif  // SRC_IGA_BASIS_FUNCTION_HANDLER_H_
