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

class AnIntegrationPoint : public Test {
 public:
  AnIntegrationPoint() : integration_point_(1.5, 0.5) {}

 protected:
  iga::itg::IntegrationPoint integration_point_;
};

TEST_F(AnIntegrationPoint, ReturnsCorrectCoordinate) { // NOLINT
  ASSERT_THAT(integration_point_.GetCoordinate(), DoubleEq(1.5));
}

TEST_F(AnIntegrationPoint, ReturnsCorrectWeight) { // NOLINT
  ASSERT_THAT(integration_point_.GetWeight(), DoubleEq(0.5));
}
