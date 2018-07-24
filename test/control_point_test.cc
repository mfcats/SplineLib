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

#include "control_point.h"

using testing::Test;
using testing::DoubleEq;

class AControlPoint : public Test {
 public:
  AControlPoint() : control_point({1.0, 2.0}) {}

 protected:
  baf::ControlPoint control_point;
};

TEST_F(AControlPoint, ReturnsCorrectDimension) { // NOLINT
  ASSERT_THAT(control_point.GetDimension(), 2);
}

TEST_F(AControlPoint, Returns1For0Dimension) { // NOLINT
  ASSERT_THAT(control_point.GetValue(0), DoubleEq(1.0));
}

TEST_F(AControlPoint, Returns2For1Dimension) { // NOLINT
  ASSERT_THAT(control_point.GetValue(1), DoubleEq(2.0));
}
