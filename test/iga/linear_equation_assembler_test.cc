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

#include <armadillo>

#include "element_integral_calculator.h"
#include "linear_equation_assembler.h"
#include "matlab_test_data_2.h"
#include "matrix.h"
#include "test_spline.h"
#include "two_point_gauss_legendre.h"

using testing::DoubleNear;

TEST_F(AnIGATestSpline, TestLeftSide) { // NOLINT
  iga::LinearEquationAssembler linear_equation_assembler = iga::LinearEquationAssembler(nurbs_);
  iga::ElementIntegralCalculator elm_itg_calc = iga::ElementIntegralCalculator(nurbs_);
  int n = nurbs_->GetNumberOfControlPoints();
  std::shared_ptr<arma::dmat> matA = std::make_shared<arma::dmat>(n, n, arma::fill::zeros);
  iga::itg::IntegrationRule rule = iga::itg::TwoPointGaussLegendre();
  linear_equation_assembler.GetLeftSide(rule, matA, elm_itg_calc);
  for (uint64_t i = 0; i < matlab_matrix_a.size(); ++i) {
    for (uint64_t j = 0; j < matlab_matrix_a[0].size(); ++j) {
      ASSERT_THAT((*matA)(i, j), DoubleNear(matlab_matrix_a[i][j], 0.0005));
    }
  }
}

TEST_F(AnIGATestSpline, TestRightSide) { // NOLINT
  iga::LinearEquationAssembler linear_equation_assembler = iga::LinearEquationAssembler(nurbs_);
  iga::ElementIntegralCalculator elm_itg_calc = iga::ElementIntegralCalculator(nurbs_);
  int n = nurbs_->GetNumberOfControlPoints();
  std::shared_ptr<arma::dvec> vecB = std::make_shared<arma::dvec>(n, arma::fill::zeros);
  std::shared_ptr<arma::dvec> srcCp = std::make_shared<arma::dvec>(n, arma::fill::ones);
  iga::itg::IntegrationRule rule = iga::itg::TwoPointGaussLegendre();
  linear_equation_assembler.GetRightSide(rule, vecB, elm_itg_calc, srcCp);
  for (uint64_t i = 0; i < matlab_vector_b.size(); ++i) {
      ASSERT_THAT((*vecB)(i), DoubleNear(matlab_vector_b[i], 0.0005));
  }
}

TEST_F(AnIGATestSpline, TestEquationSystemWithBC) { // NOLINT
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
  for (uint64_t i = 0; i < matlab_matrix_a_bc.size(); ++i) {
    for (uint64_t j = 0; j < matlab_matrix_a_bc[0].size(); ++j) {
      ASSERT_THAT((*matA)(i, j), DoubleNear(matlab_matrix_a_bc[i][j], 0.0005));
    }
  }
  for (uint64_t i = 0; i < matlab_vector_b_bc.size(); ++i) {
    ASSERT_THAT((*vecB)(i), DoubleNear(matlab_vector_b_bc[i], 0.0005));
  }
}
