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

#include "knot_vector.h"

using testing::Test;
using testing::Eq;
using testing::DoubleEq;

// A knot vector is defined as a sequence of non-decreasing real numbers (knots).
class AKnotVector : public Test {
 public:
  AKnotVector() : knot_vector_({ParamCoord{0.0}, ParamCoord{0.0}, ParamCoord{0.0}, ParamCoord{0.5}, ParamCoord{0.5}, ParamCoord{0.75}, ParamCoord{1.0}, ParamCoord{1.0}, ParamCoord{1.0}}) {}

 protected:
  baf::KnotVector knot_vector_;
};

// The i-th knot span is defined as the half-open interval (u_i, u_{i+1}]. For
// U = {0.0, 0.0, 0.0, 0.5, 0.5, 0.75, 1.0, 1.0, 1.0} u = 0.0 is in knot span 2.
// The parametric coordinate is equal to the smallest knot
TEST_F(AKnotVector, FindsZeroInKnotSpanTwo) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParamCoord{0.0}), Eq(2));
}

// The parametric coordinate is between two knots.
TEST_F(AKnotVector, Finds0_3InKnotSpanTwo) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParamCoord{0.3}), Eq(2));
}

// The parametric coordinate is last knots.
TEST_F(AKnotVector, Finds1_0InKnotSpanFive) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParamCoord{1.0}), Eq(5));
}

// The parametric coordinate is equal to an inner repeated knot.
TEST_F(AKnotVector, Finds0_5InKnotSpanFour) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParamCoord{0.5}), Eq(4));
}

// The parametric coordinate is equal to an inner not repeated knot.
TEST_F(AKnotVector, Finds0_75InKnotSpanFive) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParamCoord{0.75}), Eq(5));
}

// The last knot is defined to be in the last non-zero knot span.
TEST_F(AKnotVector, FindsLastKnotInSpanSix) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(knot_vector_.GetLastKnot()), Eq(5));
}

TEST_F(AKnotVector, CanCheckIfParametricCoordinateIsEqualLastKnot) { // NOLINT
  ASSERT_THAT(knot_vector_.IsLastKnot(knot_vector_.GetLastKnot()), Eq(true));
}

TEST_F(AKnotVector, CanCheckIfCoordinateIsNotEqualLastKnot) { // NOLINT
  ASSERT_THAT(knot_vector_.IsLastKnot(ParamCoord{0.9}), Eq(false));
}

TEST_F(AKnotVector, ReturnsCorrectKnot) { // NOLINT
  ASSERT_THAT(knot_vector_.knot(5).get(), DoubleEq(0.75));
}

TEST_F(AKnotVector, CanBeChangedWithAccessOperator) { // NOLINT
  knot_vector_[0] = ParamCoord{7.0};
  ASSERT_THAT(knot_vector_[0].get(), DoubleEq(7.0));
}

TEST_F(AKnotVector, FindsParametricCoordinateInKnotVectorRange) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInKnotVectorRange(ParamCoord{0.4}), Eq(true));
}

TEST_F(AKnotVector, FindsSmallestKnotInKnotVectorRange) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInKnotVectorRange(ParamCoord{0.0}), Eq(true));
}

TEST_F(AKnotVector, FindsLargestKnotInKnotVectorRange) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInKnotVectorRange(ParamCoord{1.0}), Eq(true));
}

TEST_F(AKnotVector, DoesNotFindSmallParametricCoordinateInKnotVectorRange) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInKnotVectorRange(ParamCoord{-0.4}), Eq(false));
}

TEST_F(AKnotVector, DoesNotFindLargeParametricCoordinateInKnotVectorRange) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInKnotVectorRange(ParamCoord{1.5}), Eq(false));
}

TEST_F(AKnotVector, CanBeCopied) { // NOLINT
  baf::KnotVector knotVector = this->knot_vector_;
  ASSERT_THAT(knotVector, Eq(this->knot_vector_));
}

TEST_F(AKnotVector, CanBeAssigned) { // NOLINT
  baf::KnotVector knotVector;
  knotVector = this->knot_vector_;
  ASSERT_THAT(knotVector, Eq(this->knot_vector_));
}

TEST_F(AKnotVector, CanBeMovedInAssignment) { // NOLINT
  baf::KnotVector knotVector;
  knotVector = baf::KnotVector({ParamCoord{0.0}, ParamCoord{0.25}, ParamCoord{0.5}});
  ASSERT_THAT(knotVector, Eq(baf::KnotVector({ParamCoord{0.0}, ParamCoord{0.25}, ParamCoord{0.5}})));
}

TEST_F(AKnotVector, CanBeAssignedByIterators) { // NOLINT
  baf::KnotVector knotVector = baf::KnotVector(knot_vector_.begin(), knot_vector_.end());
  ASSERT_THAT(knotVector, Eq(this->knot_vector_));
}
