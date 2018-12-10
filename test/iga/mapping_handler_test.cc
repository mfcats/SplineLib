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

#include "mapping_handler.h"
#include "test_spline.h"

using testing::DoubleNear;

TEST_F(AnIGATestSpline, TestMappingHandler) { // NOLINT
  iga::MappingHandler<2> mapping_handler(nurbs_);
  double j = mapping_handler.GetJacobianDeterminant(std::array<ParamCoord, 2>({ParamCoord{0.367}, ParamCoord{0.893}}));
  ASSERT_THAT(j, DoubleNear(0.0905, 0.00005));
}
