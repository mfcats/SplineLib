/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "linear_equation_assembler.h"

iga::LinearEquationAssembler::LinearEquationAssembler(std::shared_ptr<spl::NURBS<2>> spl) : spline_(std::move(spl)) {
  elm_gen_ = std::make_shared<iga::elm::ElementGenerator<2>>(spline_);
}

void iga::LinearEquationAssembler::GetLeftSide(const iga::itg::IntegrationRule &rule,
                                               const std::shared_ptr<arma::dmat> &matA,
                                               const iga::ElementIntegralCalculator &elm_itg_calc) const {
  int num_elements = static_cast<int>(elm_gen_->GetElementList(0).size() * elm_gen_->GetElementList(1).size());
  for (int e = 0; e < num_elements; ++e) {
    elm_itg_calc.GetLaplaceElementIntegral(e, rule, matA);
  }
}

void iga::LinearEquationAssembler::GetRightSide(const iga::itg::IntegrationRule &rule,
                                                const std::shared_ptr<arma::dvec> &vecB,
                                                const iga::ElementIntegralCalculator &elm_itg_calc,
                                                const std::shared_ptr<arma::dvec> &srcCp) const {
  int num_elements = static_cast<int>(elm_gen_->GetElementList(0).size() * elm_gen_->GetElementList(1).size());
  for (int e = 0; e < num_elements; ++e) {
    elm_itg_calc.GetLaplaceElementIntegral(e, rule, vecB, srcCp);
  }
}

void iga::LinearEquationAssembler::SetZeroBC(const std::shared_ptr<arma::dmat> &matA,
                                             const std::shared_ptr<arma::dvec> &vecB) {
  uint64_t l = 0;
  int n = spline_->GetPointsPerDirection()[0];
  int m = spline_->GetPointsPerDirection()[1];
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      if (!((i > 0) && (j > 0) && (i < n - 1) && (j < m - 1))) {
        (*vecB)(l) = 0;
        (*matA).row(l).fill(0);
        (*matA)(l, l) = 1;
      }
      l += 1;
    }
  }
}

// Only used for test case in test/solution_vtk_writer_examples.cc which is currently commented out.
/*void iga::LinearEquationAssembler::SetLinearBC(const std::shared_ptr<arma::dmat> &matA,
                                               const std::shared_ptr<arma::dvec> &vecB) {
  uint64_t l = 0;
  uint64_t k = 0;
  int n = spline_->GetPointsPerDirection()[0];
  int m = spline_->GetPointsPerDirection()[1];
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      if (i == 0) {
        (*vecB)(l) = 0;
        (*matA).row(l).fill(0);
        (*matA)(l, l) = 1;
      } else if (i == n - 1) {
        (*vecB)(l) = 1;
        (*matA).row(l).fill(0);
        (*matA)(l, l) = 1;
      } else if ((j == 0) || (j == m - 1)) {
        (*vecB)(l) = 0 + ((1.0 - 0.0) / (2.0 - 0.0)) * (spline_->GetControlPoints()[k] - 0.0);
        (*matA).row(l).fill(0);
        (*matA)(l, l) = 1;
      }
      l += 1;
      k += 3;
    }
  }
}*/
