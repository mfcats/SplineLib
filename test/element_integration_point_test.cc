/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "element_integration_point.h"
#include "gmock/gmock.h"
#include <vector>

using testing::DoubleEq;
using testing::Test;

class AElementIntegrationPoint : public Test{
 public:
  AElementIntegrationPoint() : basis_functions_({2.3, 4.5}), element_integration_point(basis_functions_) {}

 protected:
  std::vector<double> basis_functions_;
  elm::ElementIntegrationPoint element_integration_point;
};

TEST_F(AElementIntegrationPoint, getNonZeroBasisFunctions) { //NOLINT
  ASSERT_THAT(element_integration_point.GetNonZeroBasisFunctions()[0] , DoubleEq(2.3));
  ASSERT_THAT(element_integration_point.GetNonZeroBasisFunctions()[1] , DoubleEq(4.5));
}

TEST_F(AElementIntegrationPoint, size2) { //NOLINT
  ASSERT_THAT(element_integration_point.GetNumberOfNonZeroBasisFunctions(), 2);
}

TEST_F(AElementIntegrationPoint, valueAt0) { //NOLINT
  ASSERT_THAT(element_integration_point.GetBasisFunctionValue(1), DoubleEq(4.5));
}
