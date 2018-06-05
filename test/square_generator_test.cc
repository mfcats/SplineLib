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

#include "numeric_settings.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class ASquare : public Test {
 public:
  ASquare() {
    spl::SquareGenerator squareGenerator = spl::SquareGenerator();
    square_ = squareGenerator.CreateSquare();
  }

 protected:
  std::unique_ptr<spl::BSpline<2>> square_;
};

TEST_F(ASquare, ReturnsCorrectDegree) {
  ASSERT_THAT(square_->GetDegree(0), 2);
  ASSERT_THAT(square_->GetDegree(1), 2);
}

TEST_F(ASquare, ReturnsCorrectKnotVectorSizes) {
  ASSERT_THAT(square_->GetKnotVector(0).NumberOfKnots(), 6);
  ASSERT_THAT(square_->GetKnotVector(1).NumberOfKnots(), 6);
}

TEST_F(ASquare, ReturnsCorrectLeftLowerCorner) {
  ASSERT_THAT(square_->Evaluate({ParamCoord{0}, ParamCoord{0}}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(square_->Evaluate({0, 0}, {1})[0], DoubleEq(-1));
}

TEST_F(ASquare, ReturnsCorrectLeftUpperCorner) {
  ASSERT_THAT(square_->Evaluate({0, 1}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(square_->Evaluate({0, 1}, {1})[0], DoubleEq(1));
}

TEST_F(ASquare, ReturnsCorrectMiddle) {
  ASSERT_THAT(square_->Evaluate({0.5, 0.5}, {0})[0], DoubleEq(0));
  ASSERT_THAT(square_->Evaluate({0.5, 0.5}, {1})[0], DoubleEq(0));
}

TEST_F(ASquare, ReturnsCorrectCoordinatesOfRandomPoint) {
  ASSERT_THAT(square_->Evaluate({0.2, 0.75}, {0})[0], DoubleEq(-0.6));
  ASSERT_THAT(square_->Evaluate({0.2, 0.75}, {1})[0], DoubleEq(0.5));
}

TEST_F(ASquare, ReturnsCorrectRightLowerCorner) {
  ASSERT_THAT(square_->Evaluate({1, 0}, {0})[0], DoubleEq(1));
  ASSERT_THAT(square_->Evaluate({1, 0}, {1})[0], DoubleEq(-1));
}

TEST_F(ASquare, ReturnsCorrectRightUpperCorner) {
  ASSERT_THAT(square_->Evaluate({1, 1}, {0})[0], DoubleEq(1));
  ASSERT_THAT(square_->Evaluate({1, 1}, {1})[0], DoubleEq(1));
}

class ASquareWithDegree3And8Knots : public Test {
 public:
  ASquareWithDegree3And8Knots() {
    spl::SquareGenerator squareGenerator = spl::SquareGenerator(3, 8);
    square_ = squareGenerator.CreateSquare();
  }

 protected:
  std::unique_ptr<spl::BSpline<2>> square_;
};

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectDegree) {
  ASSERT_THAT(square_->GetDegree(0), 3);
  ASSERT_THAT(square_->GetDegree(1), 3);
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectKnotVectorSizes) {
  ASSERT_THAT(square_->GetKnotVector(0).NumberOfKnots(), 8);
  ASSERT_THAT(square_->GetKnotVector(1).NumberOfKnots(), 8);
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectLeftLowerCorner) {
  ASSERT_THAT(square_->Evaluate({0, 0}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(square_->Evaluate({0, 0}, {1})[0], DoubleEq(-1));
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectLeftUpperCorner) {
  ASSERT_THAT(square_->Evaluate({0, 1}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(square_->Evaluate({0, 1}, {1})[0], DoubleEq(1));
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectMiddle) {
  ASSERT_THAT(square_->Evaluate({0.5, 0.5}, {0})[0], DoubleNear(0, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(square_->Evaluate({0.5, 0.5}, {1})[0], DoubleNear(0, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectCoordinatesOfRandomPoint) {
  ASSERT_THAT(square_->Evaluate({0.2, 0.75}, {0})[0], DoubleEq(-0.6));
  ASSERT_THAT(square_->Evaluate({0.2, 0.75}, {1})[0], DoubleEq(0.5));
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectRightLowerCorner) {
  ASSERT_THAT(square_->Evaluate({1, 0}, {0})[0], DoubleEq(1));
  ASSERT_THAT(square_->Evaluate({1, 0}, {1})[0], DoubleEq(-1));
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectRightUpperCorner) {
  ASSERT_THAT(square_->Evaluate({1, 1}, {0})[0], DoubleEq(1));
  ASSERT_THAT(square_->Evaluate({1, 1}, {1})[0], DoubleEq(1));
}

class ASquareWithDegree3And10Knots : public Test {
 public:
  ASquareWithDegree3And10Knots() {
    spl::SquareGenerator squareGenerator = spl::SquareGenerator(3, 10);
    square_ = squareGenerator.CreateSquare();
  }

 protected:
  std::unique_ptr<spl::BSpline<2>> square_;
};

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectDegree) {
  ASSERT_THAT(square_->GetDegree(0), 3);
  ASSERT_THAT(square_->GetDegree(1), 3);
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectKnotVectorSizes) {
  ASSERT_THAT(square_->GetKnotVector(0).NumberOfKnots(), 10);
  ASSERT_THAT(square_->GetKnotVector(1).NumberOfKnots(), 10);
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectLeftLowerCorner) {
  ASSERT_THAT(square_->Evaluate({0, 0}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(square_->Evaluate({0, 0}, {1})[0], DoubleEq(-1));
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectLeftUpperCorner) {
  ASSERT_THAT(square_->Evaluate({0, 1}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(square_->Evaluate({0, 1}, {1})[0], DoubleEq(1));
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectMiddle) {
  ASSERT_THAT(square_->Evaluate({0.5, 0.5}, {0})[0], DoubleNear(0, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(square_->Evaluate({0.5, 0.5}, {1})[0], DoubleNear(0, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectCoordinatesOfRandomPoint) {
  ASSERT_THAT(square_->Evaluate({0.2, 0.75}, {0})[0], DoubleNear(-0.46, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(square_->Evaluate({0.2, 0.7}, {1})[0], DoubleEq(0.2845));
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectRightLowerCorner) {
  ASSERT_THAT(square_->Evaluate({1, 0}, {0})[0], DoubleEq(1));
  ASSERT_THAT(square_->Evaluate({1, 0}, {1})[0], DoubleEq(-1));
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectRightUpperCorner) {
  ASSERT_THAT(square_->Evaluate({1, 1}, {0})[0], DoubleEq(1));
  ASSERT_THAT(square_->Evaluate({1, 1}, {1})[0], DoubleEq(1));
}
