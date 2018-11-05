/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <array>

#include "gmock/gmock.h"

#include "element_integral_calculator.h"
#include "matlab_test_data.h"
#include "two_point_gauss_legendre.h"
#include "test_spline.h"

using testing::DoubleNear;

TEST_F(AnIGATestSpline, TestElementIntegralCalculator) { // NOLINT
  iga::ElementIntegralCalculator elm_itg_calc = iga::ElementIntegralCalculator(nurbs_);
  iga::itg::IntegrationRule rule = iga::itg::TwoPointGaussLegendre();
  int n = nurbs_->GetPointsPerDirection()[0] * nurbs_->GetPointsPerDirection()[1];
  std::shared_ptr<iga::Matrix> matrix = std::make_shared<iga::Matrix>(n, n);
  elm_itg_calc.GetLaplaceElementIntegral(0, rule, matrix);
  for (uint64_t i = 0; i < matlab_element_one_integral.size(); ++i) {
    for (uint64_t j = 0; j < matlab_element_one_integral[0].size(); ++j) {
      ASSERT_THAT(matrix->GetMatrixEntry(i, j), DoubleNear(matlab_element_one_integral[i][j], 0.00005));
    }
  }
}
