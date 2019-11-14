/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/baf/knot_vector.h"

using testing::Test;
using testing::Eq;
using testing::DoubleEq;

using namespace splinelib::src;

class AKnotVectorOfDegree2 : public Test {
 public:
  AKnotVectorOfDegree2()
      : knot_vector_({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}, ParametricCoordinate{0.0},
                      ParametricCoordinate{0.5}, ParametricCoordinate{0.5}, ParametricCoordinate{0.75},
                      ParametricCoordinate{1.0}, ParametricCoordinate{1.0}, ParametricCoordinate{1.0}}) {}

 protected:
  baf::KnotVector knot_vector_;
};

TEST_F(AKnotVectorOfDegree2, CanBeAssignedByIterators) { // NOLINT
  baf::KnotVector knot_vector = baf::KnotVector(knot_vector_.begin(), knot_vector_.end());
  ASSERT_THAT(knot_vector, Eq(knot_vector_));
}

TEST_F(AKnotVectorOfDegree2, CanBeCreatedWithMoveConstructor) { // NOLINT
  baf::KnotVector moved_knot_vector(std::move(knot_vector_));
  ASSERT_THAT(moved_knot_vector[5], Eq(ParametricCoordinate{0.75}));
}

TEST_F(AKnotVectorOfDegree2, CanBeCopied) { // NOLINT
  baf::KnotVector copied_knot_vector = knot_vector_;
  ASSERT_THAT(copied_knot_vector, Eq(knot_vector_));
}

TEST_F(AKnotVectorOfDegree2, CanBeAssigned) { // NOLINT
  baf::KnotVector knot_vector_to_assign{};
  knot_vector_to_assign = knot_vector_;
  ASSERT_THAT(knot_vector_to_assign, Eq(knot_vector_));
}

TEST_F(AKnotVectorOfDegree2, CanBeMoveAssigned) { // NOLINT
  baf::KnotVector knot_vector_to_assign{};
  knot_vector_to_assign = std::move(knot_vector_);
  ASSERT_THAT(knot_vector_to_assign[5], Eq(ParametricCoordinate{0.75}));
}

TEST_F(AKnotVectorOfDegree2, ReturnsKnot0_75AtIndex5WithMethodGetKnot) { // NOLINT
  ASSERT_THAT(knot_vector_[5], Eq(ParametricCoordinate{0.75}));
}

TEST_F(AKnotVectorOfDegree2, ReturnsKnot0_75AtIndex5WithOperator) { // NOLINT
  ASSERT_THAT(knot_vector_[5], Eq(ParametricCoordinate{0.75}));
}

TEST_F(AKnotVectorOfDegree2, SumOfAllKnotsOf4_75CanBeComputedUsingConstRangeBasedForLoop) { // NOLINT
  double sum_of_all_knots = 0;
  for (const auto &knot : knot_vector_) {
    sum_of_all_knots += knot.Get();
  }
  ASSERT_THAT(sum_of_all_knots, DoubleEq(4.75));
}

TEST_F(AKnotVectorOfDegree2, // NOLINT
       AllKnotValuesCanBeChangedTo4_0AndTheSumOfAllKnotsOf36CanBeComputedUsingNonConstRangeBasedForLoop) {
  for (auto &knot : knot_vector_) {
    knot = ParametricCoordinate{4.0};
  }
  double sum_of_all_knots = 0;
  for (const auto &knot : knot_vector_) {
    sum_of_all_knots += knot.Get();
  }
  ASSERT_THAT(sum_of_all_knots, DoubleEq(36));
}

TEST_F(AKnotVectorOfDegree2, KnotSpanOfFirstKnot0_0IsTwo) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParametricCoordinate{0.0}), Eq(KnotSpan{2}));
}

TEST_F(AKnotVectorOfDegree2, KnotSpanOf0_3BetwennTwoKnotsIsTwo) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParametricCoordinate{0.3}), Eq(KnotSpan{2}));
}

TEST_F(AKnotVectorOfDegree2, KnotSpanOfLastKnot1_0IsFive) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParametricCoordinate{1.0}), Eq(KnotSpan{5}));
}

TEST_F(AKnotVectorOfDegree2, KnotSpanOfInnerRepeatedKnot0_5IsFour) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParametricCoordinate{0.5}), Eq(KnotSpan{4}));
}

TEST_F(AKnotVectorOfDegree2, KnotSpanOf0_75IsFive) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParametricCoordinate{0.75}), Eq(KnotSpan{5}));
}

TEST_F(AKnotVectorOfDegree2, IsLastKnotReturnsTrueForParametricCoordinate1_0) { // NOLINT
  ASSERT_THAT(knot_vector_.IsLastKnot(ParametricCoordinate{1.0}), Eq(true));
}

TEST_F(AKnotVectorOfDegree2, IsLastKnotReturnsFalseForParametricCoordinate0_9) { // NOLINT
  ASSERT_THAT(knot_vector_.IsLastKnot(ParametricCoordinate{0.9}), Eq(false));
}

TEST_F(AKnotVectorOfDegree2, ParametricCoordinate0_4IsInKnotVectorRange) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInRange(ParametricCoordinate{0.4}), Eq(true));
}

TEST_F(AKnotVectorOfDegree2, SmallestParametricCoordinate0_0IsInKnotVectorRange) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInRange(ParametricCoordinate{0.0}), Eq(true));
}

TEST_F(AKnotVectorOfDegree2, LargestParametricCoordinate1_0IsInKnotVectorRange) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInRange(ParametricCoordinate{1.0}), Eq(true));
}

TEST_F(AKnotVectorOfDegree2, Minus0_4IsNotInKnotVectorRange0_0To1_0) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInRange(ParametricCoordinate{-0.4}), Eq(false));
}

TEST_F(AKnotVectorOfDegree2, 1_5IsNotInKnotVectorRange0_0To1_0) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInRange(ParametricCoordinate{1.5}), Eq(false));
}

TEST_F(AKnotVectorOfDegree2, InsertionOfKnot0_5IncreasesNumberOfKnots) {  // NOLINT
  baf::KnotVector knot_vector_copy = baf::KnotVector(knot_vector_.begin(), knot_vector_.end());
  knot_vector_copy.InsertKnot(ParametricCoordinate{0.5});
  ASSERT_THAT(knot_vector_copy.GetNumberOfKnots(), knot_vector_.GetNumberOfKnots() + 1);
  ASSERT_THAT(knot_vector_copy.GetKnotSpan(ParametricCoordinate{0.5}).Get(),
              knot_vector_.GetKnotSpan(ParametricCoordinate{0.5}).Get() + 1);
}

TEST_F(AKnotVectorOfDegree2, RemovalOfKnot0_5DecreasesNumberOfKnots) {  // NOLINT
  baf::KnotVector knot_vector_copy = baf::KnotVector(knot_vector_.begin(), knot_vector_.end());
  knot_vector_copy.RemoveKnot(ParametricCoordinate{0.5});
  ASSERT_THAT(knot_vector_copy.GetNumberOfKnots(), knot_vector_.GetNumberOfKnots() - 1);
  ASSERT_THAT(knot_vector_copy.GetKnotSpan(ParametricCoordinate{0.5}).Get(),
              knot_vector_.GetKnotSpan(ParametricCoordinate{0.5}).Get() - 1);
}
