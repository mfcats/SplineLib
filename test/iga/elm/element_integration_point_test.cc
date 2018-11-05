/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <vector>

#include "gmock/gmock.h"

#include "element_integration_point.h"

using testing::DoubleEq;
using testing::Test;

class AnElementIntegrationPoint : public Test{
 public:
  AnElementIntegrationPoint() : element_integration_point({2.3, 4.5}, 0.75) {}

 protected:
  iga::elm::ElementIntegrationPoint element_integration_point;
};

TEST_F(AnElementIntegrationPoint, GetNonZeroBasisFunctions) { //NOLINT
  ASSERT_THAT(element_integration_point.GetNonZeroBasisFunctions()[0] , DoubleEq(2.3));
  ASSERT_THAT(element_integration_point.GetNonZeroBasisFunctions()[1] , DoubleEq(4.5));
}

TEST_F(AnElementIntegrationPoint, SizeEqualsTwo) { //NOLINT
  ASSERT_THAT(element_integration_point.GetNumberOfNonZeroBasisFunctions(), 2);
}

TEST_F(AnElementIntegrationPoint, GetValueAtZero) { //NOLINT
  ASSERT_THAT(element_integration_point.GetBasisFunctionValue(1), DoubleEq(4.5));
}

TEST_F(AnElementIntegrationPoint, GetWeight) { //NOLINT
  ASSERT_THAT(element_integration_point.GetWeight(), DoubleEq(0.75));
}
