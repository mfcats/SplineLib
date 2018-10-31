/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_ELEMENT_INTEGRAL_CALCULATOR_H_
#define SRC_ELEMENT_INTEGRAL_CALCULATOR_H_

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
    mapping_handler_ = std::make_shared<iga::MappingHandler>(spline_);
    iga::ConnectivityHandler connectivity_handler(spline_);
    connectivity_ = connectivity_handler.GetConnectivity();
  }

  double GetLaplaceElementIntegral(int element_number, const iga::itg::IntegrationRule<2> &rule, iga::Matrix matA) {
    std::array<std::vector<iga::elm::ElementIntegrationPoint>, 2> elm_intgr_pnts =
        baf_handler_->EvaluateDrDxAtEveryElemIntgPnt(element_number, rule);

    std::vector<iga::itg::IntegrationPoint> intg_pnts = rule.GetIntegrationPoints();

    for (auto &intg_pnt : elm_intgr_pnts[0]) {
      for (int j = intg_pnt.GetNumberOfNonZeroBasisFunctions() - 1; j >= 0; ++j) {
        for (int k = intg_pnt.GetNumberOfNonZeroBasisFunctions() - 1; k >= 0; ++k) {

          double temp = intg_pnt * (elm_intgr_pnts[0][j]. * elm_intgr_pnts[0][k]
              + elm_intgr_pnts[1][j] * elm_intgr_pnts[1][k]) * jac_det;
          matA.AddToMatrixEntry(connectivity_(element_number, j), connectivity_(element_number, k), temp);
        }
      }
    }
  }

 private:
  std::shared_ptr<spl::NURBS<2>> spline_;
  std::shared_ptr<iga::BasisFunctionHandler> baf_handler_;
  std::shared_ptr<iga::MappingHandler> mapping_handler_;
  std::vector<std::vector<double>> connectivity_;
};
}  // namespace iga

#endif  // SRC_ELEMENT_INTEGRAL_CALCULATOR_H_
