/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include <vector>

#include "gmock/gmock.h"

#include "element_integration_point.h"

using testing::DoubleEq;
using testing::Test;

class AnElementIntegrationPoint : public Test{
 public:
  AnElementIntegrationPoint() : element_integration_point({2.3, 4.5}, 0.75, 1.75) {}

 protected:
  iga::elm::ElementIntegrationPoint<1> element_integration_point;
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

TEST_F(AnElementIntegrationPoint, GetJacobianDeterminant) { //NOLINT
  ASSERT_THAT(element_integration_point.GetJacobianDeterminant(), DoubleEq(1.75));
}

class AnElementIntegrationPointWithDerivatives : public Test{
 public:
  std::array<std::vector<double>, 2> baf_ders = {std::vector<double>({2.3, 4.5}), std::vector<double>({3.2, 5.4})};
  AnElementIntegrationPointWithDerivatives() : element_integration_point(baf_ders, 0.75, 1.75) {}

 protected:
  iga::elm::ElementIntegrationPoint<2> element_integration_point;
};

TEST_F(AnElementIntegrationPointWithDerivatives, GetNonZeroBasisFunctionDerivatives) { //NOLINT
  ASSERT_THAT(element_integration_point.GetNonZeroBasisFunctionDerivatives(0)[0] , DoubleEq(2.3));
  ASSERT_THAT(element_integration_point.GetNonZeroBasisFunctionDerivatives(0)[1] , DoubleEq(4.5));
  ASSERT_THAT(element_integration_point.GetNonZeroBasisFunctionDerivatives(1)[0] , DoubleEq(3.2));
  ASSERT_THAT(element_integration_point.GetNonZeroBasisFunctionDerivatives(1)[1] , DoubleEq(5.4));
}
