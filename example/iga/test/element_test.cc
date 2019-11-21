/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include <b_spline.h>
#include "gmock/gmock.h"

#include "element.h"

using testing::Test;
using testing::DoubleEq;

class A1DElement : public Test {
 public:
  A1DElement() : element({ParametricCoordinate(0.5), ParametricCoordinate(1.0)}) {}

 protected:
  iga::elm::Element element;
};

TEST_F(A1DElement, ReturnsCorrectNode) { // NOLINT
  ASSERT_THAT(element.GetLowerBound().Get(), DoubleEq(0.5));
  ASSERT_THAT(element.GetUpperBound().Get(), DoubleEq(1.0));
}
