/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "element_integral_calculator.h"

iga::ElementIntegralCalculator::ElementIntegralCalculator(
    std::shared_ptr<spl::NURBS<2>> spl) : spline_(std::move(spl)) {
  baf_handler_ = std::make_shared<iga::BasisFunctionHandler>(spline_);
  connectivity_handler_ = std::make_shared<iga::ConnectivityHandler>(spline_);
}

void iga::ElementIntegralCalculator::GetLaplaceElementIntegral(int element_number,
                                                               const iga::itg::IntegrationRule &rule,
                                                               const std::shared_ptr<arma::dmat> &matA) const {
  std::vector<iga::elm::ElementIntegrationPoint> elm_intgr_pnts =
      baf_handler_->EvaluateAllElementNonZeroNURBSBafDerivativesPhysical(element_number, rule);
  for (auto &p : elm_intgr_pnts) {
    for (int j = 0; j < p.GetNumberOfNonZeroBasisFunctionDerivatives(1); ++j) {
      for (int k = 0; k < p.GetNumberOfNonZeroBasisFunctionDerivatives(0); ++k) {
        double temp = (p.GetBasisFunctionDerivativeValue(j, 0) * p.GetBasisFunctionDerivativeValue(k, 0)
            + p.GetBasisFunctionDerivativeValue(j, 1) * p.GetBasisFunctionDerivativeValue(k, 1))
            * p.GetWeight() * p.GetJacobianDeterminant();
        (*matA)(static_cast<uint64_t>(connectivity_handler_->GetGlobalIndex(element_number, j) - 1),
                static_cast<uint64_t>(connectivity_handler_->GetGlobalIndex(element_number, k) - 1)) += temp;
      }
    }
  }
}

void iga::ElementIntegralCalculator::GetLaplaceElementIntegral(int element_number,
                                                               const iga::itg::IntegrationRule &rule,
                                                               const std::shared_ptr<arma::dvec> &vecB,
                                                               const std::shared_ptr<arma::dvec> &srcCp) const {
  std::vector<iga::elm::ElementIntegrationPoint> elm_intgr_pnts =
      baf_handler_->EvaluateAllElementNonZeroNURBSBasisFunctions(element_number, rule);
  for (auto &p : elm_intgr_pnts) {
    double bc_int_pnt = 0;
    for (int j = 0; j < p.GetNumberOfNonZeroBasisFunctions(); ++j) {
      bc_int_pnt += p.GetBasisFunctionValue(j) *
          (*srcCp)(static_cast<uint64_t>(connectivity_handler_->GetGlobalIndex(element_number, j) - 1));
    }
    for (int j = 0; j < p.GetNumberOfNonZeroBasisFunctions(); ++j) {
      double temp = p.GetBasisFunctionValue(j) * bc_int_pnt * p.GetWeight() * p.GetJacobianDeterminant();
      (*vecB)(static_cast<uint64_t>(connectivity_handler_->GetGlobalIndex(element_number, j) - 1)) += temp;
    }
  }
}
