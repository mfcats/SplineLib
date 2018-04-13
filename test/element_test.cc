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

#include "element.h"

using testing::Test;

class A1DElement : public Test {
 public:
  A1DElement() : element(1, {ControlPoint({1.0, 1.0}), ControlPoint({1.0, 2.0})}) {}

 protected:
  Element element;
};

TEST_F(A1DElement, ReturnsCorrectDimension) { // NOLINT
  ASSERT_THAT(element.dimension(), 1);
}

TEST_F(A1DElement, ReturnsCorrectNumberOfNodes) { // NOLINT
  ASSERT_THAT(element.numberOfNodes(), 2);
}

TEST_F(A1DElement, ReturnsCorrectNode) { // NOLINT
  ASSERT_THAT(element.node(1).dimension(), 2);
  ASSERT_THAT(element.node(1).GetValue(0), 1);
  ASSERT_THAT(element.node(1).GetValue(1), 2);
}
