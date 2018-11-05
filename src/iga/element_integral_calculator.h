/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_ELEMENT_INTEGRAL_CALCULATOR_H_
#define SRC_IGA_ELEMENT_INTEGRAL_CALCULATOR_H_

#include <vector>

#include "basis_function_handler.h"
#include "connectivity_handler.h"
#include "element_integration_point.h"
#include "integration_point.h"
#include "mapping_handler.h"
#include "matrix.h"
#include "nurbs.h"

namespace iga {
class ElementIntegralCalculator {
 public:
  explicit ElementIntegralCalculator(std::shared_ptr<spl::NURBS<2>> spl) : spline_(std::move(spl)) {
    baf_handler_ = std::make_shared<iga::BasisFunctionHandler>(spline_);
    iga::ConnectivityHandler connectivity_handler(spline_);
    connectivity_ = connectivity_handler.GetConnectivity();
  }

  void GetLaplaceElementIntegral(int element_number, const iga::itg::IntegrationRule &rule,
      const std::shared_ptr<iga::Matrix> &matA) {
    std::array<std::vector<iga::elm::ElementIntegrationPoint>, 2> elm_intgr_pnts =
        baf_handler_->EvaluateDrDxAtEveryElemItgPnt(element_number, rule);
    for (int i = 0; i < elm_intgr_pnts[0].size(); ++i) {
      for (int j = 0; j < elm_intgr_pnts[0][i].GetNumberOfNonZeroBasisFunctions(); ++j) {
        for (int k = 0; k < elm_intgr_pnts[0][i].GetNumberOfNonZeroBasisFunctions(); ++k) {
          double temp = (elm_intgr_pnts[0][i].GetBasisFunctionValue(j) * elm_intgr_pnts[0][i].GetBasisFunctionValue(k) +
              elm_intgr_pnts[1][i].GetBasisFunctionValue(j) * elm_intgr_pnts[1][i].GetBasisFunctionValue(k)) *
                  elm_intgr_pnts[0][i].GetWeight() * elm_intgr_pnts[0][i].GetJacDet();
          matA->AddToMatrixEntry(connectivity_[element_number][j] - 1, connectivity_[element_number][k] - 1, temp);
        }
      }
    }
  }

 private:
  std::shared_ptr<spl::NURBS<2>> spline_;
  std::shared_ptr<iga::BasisFunctionHandler> baf_handler_;
  std::vector<std::vector<int>> connectivity_;
};
}  // namespace iga

#endif  // SRC_IGA_ELEMENT_INTEGRAL_CALCULATOR_H_
