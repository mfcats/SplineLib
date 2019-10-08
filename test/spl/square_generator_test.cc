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

using namespace splinelib::src;

class ASquare : public Test {
 public:
  ASquare() {
    spl::SquareGenerator squareGenerator = spl::SquareGenerator();
    square_ = squareGenerator.CreateSquare();
  }

 protected:
  std::unique_ptr<spl::BSpline<2>> square_;
};

TEST_F(ASquare, ReturnsCorrectDegree) { // NOLINT
  ASSERT_THAT(square_->GetDegree(0), Degree{2});
  ASSERT_THAT(square_->GetDegree(1), Degree{2});
}

TEST_F(ASquare, ReturnsCorrectKnotVectorSizes) { // NOLINT
  ASSERT_THAT(square_->GetKnotVector(0)->GetNumberOfKnots(), 6);
  ASSERT_THAT(square_->GetKnotVector(1)->GetNumberOfKnots(), 6);
}

TEST_F(ASquare, ReturnsCorrectLeftLowerCorner) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0}, ParametricCoordinate{0}}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0}, ParametricCoordinate{0}}, {1})[0], DoubleEq(-1));
}

TEST_F(ASquare, ReturnsCorrectLeftUpperCorner) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0}, ParametricCoordinate{1}}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0}, ParametricCoordinate{1}}, {1})[0], DoubleEq(1));
}

TEST_F(ASquare, ReturnsCorrectMiddle) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0.5}, ParametricCoordinate{0.5}}, {0})[0], DoubleEq(0));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0.5}, ParametricCoordinate{0.5}}, {1})[0], DoubleEq(0));
}

TEST_F(ASquare, ReturnsCorrectCoordinatesOfRandomPoint) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0.2}, ParametricCoordinate{0.75}}, {0})[0], DoubleEq(-0.6));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0.2}, ParametricCoordinate{0.75}}, {1})[0], DoubleEq(0.5));
}

TEST_F(ASquare, ReturnsCorrectRightLowerCorner) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{1}, ParametricCoordinate{0}}, {0})[0], DoubleEq(1));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{1}, ParametricCoordinate{0}}, {1})[0], DoubleEq(-1));
}

TEST_F(ASquare, ReturnsCorrectRightUpperCorner) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{1}, ParametricCoordinate{1}}, {0})[0], DoubleEq(1));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{1}, ParametricCoordinate{1}}, {1})[0], DoubleEq(1));
}

class ASquareWithDegree3And8Knots : public Test {
 public:
  ASquareWithDegree3And8Knots() {
    spl::SquareGenerator squareGenerator = spl::SquareGenerator(Degree{3}, 8);
    square_ = squareGenerator.CreateSquare();
  }

 protected:
  std::unique_ptr<spl::BSpline<2>> square_;
};

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectDegree) { // NOLINT
  ASSERT_THAT(square_->GetDegree(0), Degree{3});
  ASSERT_THAT(square_->GetDegree(1), Degree{3});
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectKnotVectorSizes) { // NOLINT
  ASSERT_THAT(square_->GetKnotVector(0)->GetNumberOfKnots(), 8);
  ASSERT_THAT(square_->GetKnotVector(1)->GetNumberOfKnots(), 8);
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectLeftLowerCorner) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0}, ParametricCoordinate{0}}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0}, ParametricCoordinate{0}}, {1})[0], DoubleEq(-1));
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectLeftUpperCorner) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0}, ParametricCoordinate{1}}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0}, ParametricCoordinate{1}}, {1})[0], DoubleEq(1));
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectMiddle) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0.5}, ParametricCoordinate{0.5}}, {0})[0],
              DoubleNear(0, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0.5}, ParametricCoordinate{0.5}}, {1})[0],
              DoubleNear(0, util::numeric_settings::GetEpsilon<double>()));
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectCoordinatesOfRandomPoint) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0.2}, ParametricCoordinate{0.75}}, {0})[0], DoubleEq(-0.6));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0.2}, ParametricCoordinate{0.75}}, {1})[0], DoubleEq(0.5));
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectRightLowerCorner) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{1}, ParametricCoordinate{0}}, {0})[0], DoubleEq(1));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{1}, ParametricCoordinate{0}}, {1})[0], DoubleEq(-1));
}

TEST_F(ASquareWithDegree3And8Knots, ReturnsCorrectRightUpperCorner) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{1}, ParametricCoordinate{1}}, {0})[0], DoubleEq(1));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{1}, ParametricCoordinate{1}}, {1})[0], DoubleEq(1));
}

class ASquareWithDegree3And10Knots : public Test {
 public:
  ASquareWithDegree3And10Knots() {
    spl::SquareGenerator squareGenerator = spl::SquareGenerator(Degree{3}, 10);
    square_ = squareGenerator.CreateSquare();
  }

 protected:
  std::unique_ptr<spl::BSpline<2>> square_;
};

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectDegree) { // NOLINT
  ASSERT_THAT(square_->GetDegree(0), Degree{3});
  ASSERT_THAT(square_->GetDegree(1), Degree{3});
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectKnotVectorSizes) { // NOLINT
  ASSERT_THAT(square_->GetKnotVector(0)->GetNumberOfKnots(), 10);
  ASSERT_THAT(square_->GetKnotVector(1)->GetNumberOfKnots(), 10);
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectLeftLowerCorner) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0}, ParametricCoordinate{0}}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0}, ParametricCoordinate{0}}, {1})[0], DoubleEq(-1));
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectLeftUpperCorner) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0}, ParametricCoordinate{1}}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0}, ParametricCoordinate{1}}, {1})[0], DoubleEq(1));
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectMiddle) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0.5}, ParametricCoordinate{0.5}}, {0})[0],
              DoubleNear(0, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0.5}, ParametricCoordinate{0.5}}, {1})[0],
              DoubleNear(0, util::numeric_settings::GetEpsilon<double>()));
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectCoordinatesOfRandomPoint) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0.2}, ParametricCoordinate{0.75}}, {0})[0],
              DoubleNear(-0.46, util::numeric_settings::GetEpsilon<double>()));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{0.2}, ParametricCoordinate{0.7}}, {1})[0], DoubleEq(0.2845));
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectRightLowerCorner) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{1}, ParametricCoordinate{0}}, {0})[0], DoubleEq(1));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{1}, ParametricCoordinate{0}}, {1})[0], DoubleEq(-1));
}

TEST_F(ASquareWithDegree3And10Knots, ReturnsCorrectRightUpperCorner) { // NOLINT
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{1}, ParametricCoordinate{1}}, {0})[0], DoubleEq(1));
  ASSERT_THAT(square_->Evaluate({ParametricCoordinate{1}, ParametricCoordinate{1}}, {1})[0], DoubleEq(1));
}
