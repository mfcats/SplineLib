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

#include "integration_point.h"

using testing::Test;
using testing::DoubleEq;

class A1DIntegrationPoint : public Test {
 public:
  A1DIntegrationPoint() : integration_point_({1.5}, 0.5) {}

 protected:
  IntegrationPoint<1> integration_point_;
};

TEST_F(A1DIntegrationPoint, ReturnsCorrectCoordinate) {
  ASSERT_THAT(integration_point_.GetCoordinates().size(), 1);
  ASSERT_THAT(integration_point_.GetCoordinates()[0], DoubleEq(1.5));
}

TEST_F(A1DIntegrationPoint, ReturnsCorrectWeight) {
  ASSERT_THAT(integration_point_.GetWeight(), DoubleEq(0.5));
}

TEST_F(A1DIntegrationPoint, ReturnsCorrectDimension) {
  ASSERT_THAT(integration_point_.GetDimension(), 1);
}
