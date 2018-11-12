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
#include <array>

#include "gmock/gmock.h"
#include "matlab_test_data.h"
#include "test_spline.h"

using testing::DoubleNear;

TEST_F(AnIGATestSpline, TestElementIntegralCalculator) { // NOLINT
  elm_itg_calc.GetLaplaceElementIntegral(0, rule, matA);
  for (uint64_t i = 0; i < matlab_element_one_integral.size(); ++i) {
    for (uint64_t j = 0; j < matlab_element_one_integral[0].size(); ++j) {
      ASSERT_THAT((*matA)(i, j), DoubleNear(matlab_element_one_integral[i][j], 0.00005));
    }
  }
}
