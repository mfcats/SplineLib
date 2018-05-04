/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "square_generator.h"

#include "gmock/gmock.h"

using testing::Test;
using testing::DoubleEq;

class ASquare : public Test {
 public:
  ASquare() {
    SquareGenerator squareGenerator = SquareGenerator();
    square_ = squareGenerator.CreateSquare();
  }

 protected:
  std::unique_ptr<BSpline<2>> square_;
};

TEST_F(ASquare, ReturnsCorrectDegree) {
  ASSERT_THAT(square_->GetDegree(0), 2);
  ASSERT_THAT(square_->GetDegree(1), 2);
}

TEST_F(ASquare, ReturnsCorrectKnotVectorSizes) {
  ASSERT_THAT(square_->GetKnotVector(0).Size(), 6);
  ASSERT_THAT(square_->GetKnotVector(1).Size(), 6);
}

TEST_F(ASquare, ReturnsCorrectEdge) {
  ASSERT_THAT(square_->Evaluate({0, 0}, {0})[0], DoubleEq(-1));
}

TEST_F(ASquare, ReturnsCorrectMiddle) {
  ASSERT_THAT(square_->Evaluate({0.5, 0.5}, {0})[0], DoubleEq(0));
}

TEST_F(ASquare, ReturnsCorrectLastEdge) {
  ASSERT_THAT(square_->Evaluate({1, 1}, {0})[0], DoubleEq(1));
}

class ASquareWithDegree3And8Knots : public Test {
 public:
  ASquareWithDegree3And8Knots() {
    SquareGenerator squareGenerator = SquareGenerator(3, 8);
    square_ = squareGenerator.CreateSquare();
  }

 protected:
  std::unique_ptr<BSpline<2>> square_;
};

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectDegree) {
  ASSERT_THAT(square_->GetDegree(0), 3);
  ASSERT_THAT(square_->GetDegree(1), 3);
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectKnotVectorSizes) {
  ASSERT_THAT(square_->GetKnotVector(0).Size(), 8);
  ASSERT_THAT(square_->GetKnotVector(1).Size(), 8);
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectControlPoints) {
}

class ASquareWithDegree3And10Knots : public Test {
 public:
  ASquareWithDegree3And10Knots() {
    SquareGenerator squareGenerator = SquareGenerator(3, 10);
    square_ = squareGenerator.CreateSquare();
  }

 protected:
  std::unique_ptr<BSpline<2>> square_;
};

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectDegree) {
  ASSERT_THAT(square_->GetDegree(0), 3);
  ASSERT_THAT(square_->GetDegree(1), 3);
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectKnotVectorSizes) {
  ASSERT_THAT(square_->GetKnotVector(0).Size(), 10);
  ASSERT_THAT(square_->GetKnotVector(1).Size(), 10);
}