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

#include <array>
#include <vector>

#include "element_integration_point.h"
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
        double temp = 0;
        temp = basis_functions[0][j] * basis_functions[1][i] * spline_->GetWeight({j, i});
        sum += temp;
        nurbs_basis_functions.emplace_back(temp);
      }
    }
    for (int i = 0; i < nurbs_basis_functions.size(); ++i) {
      nurbs_basis_functions[i] = nurbs_basis_functions[i] / sum;
    }
    return nurbs_basis_functions;
  }

 private:
  std::shared_ptr<spl::NURBS<2>> spline_;
};
}  // namespace iga

#endif  // SRC_IGA_BASIS_FUNCTION_HANDLER_H_
