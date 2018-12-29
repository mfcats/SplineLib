/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_BDF_HANDLER_H_
#define SRC_IGA_BDF_HANDLER_H_

#include "linear_equation_assembler.h"

namespace iga {
template<int DIM>
class BDFHandler {
 public:
  explicit BDFHandler(std::shared_ptr<spl::NURBS<DIM>> spl) : spline_(std::move(spl)) {
    elm_gen_ = std::make_shared<iga::elm::ElementGenerator<DIM>>(spline_);
    baf_handler_ = std::make_shared<iga::BasisFunctionHandler<DIM>>(spline_);
    connectivity_handler_ = std::make_shared<iga::ConnectivityHandler<DIM>>(spline_);
  }

  std::shared_ptr<arma::dmat> GetBDF1LeftSide(const iga::itg::IntegrationRule &rule,
                                              const std::shared_ptr<arma::dmat> &matA, double dt) {
    auto left = std::make_shared<arma::dmat>((*matA) + (GetTimeDiscretizationMatrix(rule) / dt));
    return left;
  }

  std::shared_ptr<arma::dvec> GetBDF1RightSide(const iga::itg::IntegrationRule &rule,
                                               const std::shared_ptr<arma::dvec> &vecB,
                                               const std::shared_ptr<arma::dvec> &prevSol, double dt) {
    auto right = std::make_shared<arma::dvec>((*vecB) + (GetTimeDiscretizationMatrix(rule) / dt) * (*prevSol));
    return right;
  }

 private:
  arma::dmat GetTimeDiscretizationMatrix(const iga::itg::IntegrationRule &rule) const {
    auto num_cp = static_cast<uint64_t>(spline_->GetNumberOfControlPoints());
    arma::dmat matA2(num_cp, num_cp, arma::fill::zeros);
    for (int e = 0; e < elm_gen_->GetNumberOfElements(); ++e) {
      std::vector<iga::elm::ElementIntegrationPoint<DIM>> elm_intgr_pnts =
          baf_handler_->EvaluateAllElementNonZeroNURBSBasisFunctions(e, rule);
      for (auto &p : elm_intgr_pnts) {
        for (int j = 0; j < p.GetNumberOfNonZeroBasisFunctions(); ++j) {
          for (int k = 0; k < p.GetNumberOfNonZeroBasisFunctions(); ++k) {
            double temp = p.GetBasisFunctionValue(j) * p.GetBasisFunctionValue(k) * p.GetWeight()
                * p.GetJacobianDeterminant();
            matA2(static_cast<uint64_t>(connectivity_handler_->GetGlobalIndex(e, j) - 1),
                    static_cast<uint64_t>(connectivity_handler_->GetGlobalIndex(e, k) - 1)) += temp;
          }
        }
      }
    }
    return matA2;
  }

  std::shared_ptr<spl::NURBS<DIM>> spline_;
  std::shared_ptr<iga::elm::ElementGenerator<DIM>> elm_gen_;
  std::shared_ptr<iga::ConnectivityHandler<DIM>> connectivity_handler_;
  std::shared_ptr<iga::BasisFunctionHandler<DIM>> baf_handler_;
};
}

#endif  // SRC_IGA_BDF_HANDLER_H_
