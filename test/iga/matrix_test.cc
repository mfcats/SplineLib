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

#include "matrix.h"

using testing::Test;
using testing::DoubleEq;

class AMatrix : public Test {
 public:
  AMatrix() {}
};

TEST_F(AMatrix, TestMatrix) { // NOLINT
  iga::Matrix matA(2, 2);
  matA.WriteToMatrix(1, 1, 2.5);
  ASSERT_THAT(matA.GetMatrixEntry(1, 1), DoubleEq(2.5));
  ASSERT_THAT(matA.GetMatrixEntry(0, 0), DoubleEq(0));
  ASSERT_THAT(matA.GetMatrixEntry(0, 1), DoubleEq(0));
  ASSERT_THAT(matA.GetMatrixEntry(1, 0), DoubleEq(0));
  matA.AddToMatrixEntry(0, 1, 1.5);
  ASSERT_THAT(matA.GetMatrixEntry(1, 1), DoubleEq(2.5));
  ASSERT_THAT(matA.GetMatrixEntry(0, 0), DoubleEq(0));
  ASSERT_THAT(matA.GetMatrixEntry(0, 1), DoubleEq(1.5));
  ASSERT_THAT(matA.GetMatrixEntry(1, 0), DoubleEq(0));
}
