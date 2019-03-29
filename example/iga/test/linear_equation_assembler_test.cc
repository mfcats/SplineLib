/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <armadillo>

#include "gmock/gmock.h"
#include "matlab_test_data_2.h"
#include "test_spline.h"

using testing::DoubleNear;

TEST_F(AnIGATestSpline, TestLeftSide) { // NOLINT
  linear_equation_assembler.GetLeftSide(rule, matA, elm_itg_calc);
  for (uint64_t i = 0; i < matlab_matrix_a.size(); ++i) {
    for (uint64_t j = 0; j < matlab_matrix_a[0].size(); ++j) {
      ASSERT_THAT((*matA)(i, j), DoubleNear(matlab_matrix_a[i][j], 0.00005));
    }
  }
}

TEST_F(AnIGATestSpline, TestRightSide) { // NOLINT
  linear_equation_assembler.GetRightSide(rule, vecB, elm_itg_calc, srcCp);
  for (uint64_t i = 0; i < matlab_vector_b.size(); ++i) {
      ASSERT_THAT((*vecB)(i), DoubleNear(matlab_vector_b[i], 0.00005));
  }
}

TEST_F(AnIGATestSpline, TestEquationSystemWithBC) { // NOLINT
  linear_equation_assembler.GetLeftSide(rule, matA, elm_itg_calc);
  linear_equation_assembler.GetRightSide(rule, vecB, elm_itg_calc, srcCp);
  linear_equation_assembler.SetZeroBC(matA, vecB);
  for (uint64_t i = 0; i < matlab_matrix_a_bc.size(); ++i) {
    for (uint64_t j = 0; j < matlab_matrix_a_bc[0].size(); ++j) {
      ASSERT_THAT((*matA)(i, j), DoubleNear(matlab_matrix_a_bc[i][j], 0.00005));
    }
  }
  for (uint64_t i = 0; i < matlab_vector_b_bc.size(); ++i) {
    ASSERT_THAT((*vecB)(i), DoubleNear(matlab_vector_b_bc[i], 0.00005));
  }
}

TEST_F(AnIGATestSpline, TestSolution) { // NOLINT
  iga::itg::IntegrationRule rule = iga::itg::TwoPointGaussLegendre();
  linear_equation_assembler.GetLeftSide(rule, matA, elm_itg_calc);
  linear_equation_assembler.GetRightSide(rule, vecB, elm_itg_calc, srcCp);
  linear_equation_assembler.SetZeroBC(matA, vecB);
  arma::dvec solution = arma::solve(*matA, *vecB);
  for (uint64_t i = 0; i < matlab_solution.size(); ++i) {
    ASSERT_THAT(solution(i), DoubleNear(matlab_solution[i], 0.00005));
  }
}
