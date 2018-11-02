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

#include "linear_equation_assembler.h"
#include "matlab_test_data_2.h"
#include "matrix.h"
#include "test_spline.h"
#include "two_point_gauss_legendre.h"

using testing::DoubleNear;

TEST_F(AnIGATestSpline, TestLinearEquationAssembler) { // NOLINT
  iga::LinearEquationAssembler linear_equation_assembler = iga::LinearEquationAssembler(nurbs_);
  std::shared_ptr<iga::Matrix> matA = std::make_shared<iga::Matrix>(49, 49);
  iga::itg::IntegrationRule rule = iga::itg::TwoPointGaussLegendre();
  linear_equation_assembler.Laplace(rule, matA);
  for (int i = 0; i < matlab_matrix_a.size(); ++i) {
    for (int j = 0; j < matlab_matrix_a[0].size(); ++j) {
      ASSERT_THAT(matA->GetMatrixEntry(i, j), DoubleNear(matlab_matrix_a[i][j], 0.0005));
    }
  }
}
