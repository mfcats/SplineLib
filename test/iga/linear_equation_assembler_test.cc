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

TEST_F(AnIGATestSpline, TestLinearEquationAssembler) { // NOLINT
  iga::LinearEquationAssembler linear_equation_assembler = iga::LinearEquationAssembler(nurbs_);
  iga::ElementIntegralCalculator elm_itg_calc = iga::ElementIntegralCalculator(nurbs_);
  std::shared_ptr<arma::dmat> matA = std::make_shared<arma::dmat>(49, 49, arma::fill::zeros);
  iga::itg::IntegrationRule rule = iga::itg::TwoPointGaussLegendre();
  linear_equation_assembler.GetLeftSide(rule, matA, elm_itg_calc);
  for (uint64_t i = 0; i < matlab_matrix_a.size(); ++i) {
    for (uint64_t j = 0; j < matlab_matrix_a[0].size(); ++j) {
      ASSERT_THAT((*matA)(i, j), DoubleNear(matlab_matrix_a[i][j], 0.0005));
    }
  }
}
