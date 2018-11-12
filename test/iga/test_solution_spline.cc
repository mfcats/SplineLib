/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "gmock/gmock.h"

#include <any>
#include <armadillo>

#include "iges_writer.h"
#include "linear_equation_assembler.h"
#include "nurbs.h"
#include "solution_spline.h"
#include "test_spline.h"
#include "two_point_gauss_legendre.h"

TEST_F(AnIGATestSpline, TestSolutionSpline) { // NOLINT
  iga::LinearEquationAssembler linear_equation_assembler = iga::LinearEquationAssembler(nurbs_);
  iga::ElementIntegralCalculator elm_itg_calc = iga::ElementIntegralCalculator(nurbs_);
  int n = nurbs_->GetNumberOfControlPoints();
  std::shared_ptr<arma::dmat> matA = std::make_shared<arma::dmat>(n, n, arma::fill::zeros);
  std::shared_ptr<arma::dvec> vecB = std::make_shared<arma::dvec>(n, arma::fill::zeros);
  std::shared_ptr<arma::dvec> srcCp = std::make_shared<arma::dvec>(n, arma::fill::ones);
  iga::itg::IntegrationRule rule = iga::itg::TwoPointGaussLegendre();
  linear_equation_assembler.GetLeftSide(rule, matA, elm_itg_calc);
  linear_equation_assembler.GetRightSide(rule, vecB, elm_itg_calc, srcCp);
  linear_equation_assembler.SetZeroBC(matA, vecB);
  arma::dvec solution = arma::solve(*matA, *vecB);
  iga::SolutionSpline sol_spl(nurbs_, solution);
  std::any sol_spl_ = std::make_any<std::shared_ptr<spl::NURBS<2>>>(sol_spl.GetSolutionSpline());
  io::IGESWriter iges_writer;
  iges_writer.WriteFile({sol_spl_}, "solution_spline.iges");
  remove("solution_spline.iges");
}
