/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_POISSON_H_
#define SRC_IGA_POISSON_H_

#include <armadillo>

#include "bdf_handler.h"
#include "linear_equation_assembler.h"
#include "nurbs.h"
#include "spline.h"

namespace iga {
template<int DIM>
class PoissonProblem {
 public:
  PoissonProblem(std::shared_ptr<spl::NURBS<DIM>> spl, const iga::itg::IntegrationRule &rule) :
  spline_(std::move(spl)), num_cp_(spline_->GetNumberOfControlPoints()), rule_(rule) {
    linear_equation_assembler_ = std::make_shared<iga::LinearEquationAssembler<DIM>>(spline_);
    elm_itg_calc_ = std::make_shared<iga::ElementIntegralCalculator<DIM>>(spline_);
    matA_ = std::make_shared<arma::dmat>(num_cp_, num_cp_, arma::fill::zeros);
    vecB_ = std::make_shared<arma::dvec>(num_cp_, arma::fill::zeros);
    srcCp_ = std::make_shared<arma::dvec>(num_cp_, arma::fill::ones);
  }

  arma::dvec GetSteadyStateSolution() { // int DirichletBC, int constSrc) {
    linear_equation_assembler_->GetLeftSide(rule_, matA_, *elm_itg_calc_);
    linear_equation_assembler_->GetRightSide(rule_, vecB_, *elm_itg_calc_, srcCp_);
    linear_equation_assembler_->SetZeroBC(matA_, vecB_);
    return arma::solve(*matA_, *vecB_);
  }

  std::vector<std::shared_ptr<arma::dvec>> GetUnsteadyStateSolution(double dt, double tEnd) { // int DirichletBC, int constSrc
    iga::BDFHandler bdf_handler(spline_, rule_);
    std::vector<std::shared_ptr<arma::dvec>> solutions;
    int timeSteps = static_cast<int>(tEnd / dt);
    std::shared_ptr<arma::dvec> uprev = std::make_shared<arma::dvec>(static_cast<uint64_t>(num_cp_), arma::fill::zeros);
    solutions.emplace_back(uprev);
    linear_equation_assembler_->GetLeftSide(rule_, matA_, *elm_itg_calc_);
    linear_equation_assembler_->GetRightSide(rule_, vecB_, *elm_itg_calc_, srcCp_);
    auto left = bdf_handler.GetBDF1LeftSide(matA_, dt);
    for (int i = 1; i <= timeSteps; ++i) {
      auto right = bdf_handler.GetBDF1RightSide(vecB_, uprev, dt);
      linear_equation_assembler_->SetZeroBC(left, right);
      uprev = std::make_shared<arma::dvec>(arma::solve(*left, *right));
      solutions.emplace_back(uprev);
    }
    return solutions;
  }

 private:
  std::shared_ptr<spl::NURBS<DIM>> spline_;
  std::shared_ptr<iga::LinearEquationAssembler<DIM>> linear_equation_assembler_;
  std::shared_ptr<iga::ElementIntegralCalculator<DIM>> elm_itg_calc_;
  iga::itg::IntegrationRule rule_;
  std::shared_ptr<arma::dmat> matA_;
  std::shared_ptr<arma::dvec> vecB_;
  std::shared_ptr<arma::dvec> srcCp_;
  int num_cp_;
};
}

#endif  // SRC_IGA_POISSON_H_
