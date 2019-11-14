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
  ASSERT_THAT(knot_vector_.GetKnot(5), Eq(ParametricCoordinate{0.75}));
}

TEST_F(AKnotVectorOfDegree2, ReturnsKnot0_75AtIndex5WithOperator) { // NOLINT
  ASSERT_THAT(knot_vector_[5], Eq(ParametricCoordinate{0.75}));
}

TEST_F(AKnotVectorOfDegree2, CanBeUsedWithConstRangeBasedForLoop) { // NOLINT
  double sum_of_all_knots = 0;
  for (const auto &knot : knot_vector_) {
    sum_of_all_knots += knot.Get();
  }
  ASSERT_THAT(sum_of_all_knots, DoubleEq(4.75));
}

TEST_F(AKnotVectorOfDegree2, CanBeUsedWithNonConstRangeBasedForLoop) { // NOLINT
  for (auto &knot : knot_vector_) {
    knot = ParametricCoordinate{4.0};
  }
  double sum_of_all_knots = 0;
  for (const auto &knot : knot_vector_) {
    sum_of_all_knots += knot.Get();
  }
  ASSERT_THAT(sum_of_all_knots, DoubleEq(36));
}

// Test GetKnotSpan if parametric coordinate is equal to the smallest knot.
TEST_F(AKnotVectorOfDegree2, FindsZeroInKnotSpanTwo) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParametricCoordinate{0.0}), Eq(KnotSpan{2}))
    << "The knot span of the smallest knot has to equal degree p = 2.";
}

// Test GetKnotSpan if parametric coordinate is between two knots.
TEST_F(AKnotVectorOfDegree2, Finds0_3InKnotSpanTwo) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParametricCoordinate{0.3}), Eq(KnotSpan{2}));
}

// Test GetKnotSpan if parametric coordinate is last knot.
// The last knot is defined to be in the last non-zero knot span.
TEST_F(AKnotVectorOfDegree2, FindsLastKnot1_0InKnotSpanFive) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParametricCoordinate{1.0}), Eq(KnotSpan{5}));
}

// Test GetKnotSpan if parametric coordinate is equal to an inner repeated knot.
TEST_F(AKnotVectorOfDegree2, Finds0_5InKnotSpanFour) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParametricCoordinate{0.5}), Eq(KnotSpan{4}));
}

// Test GetKnotSpan if parametric coordinate is equal to an inner not repeated knot.
TEST_F(AKnotVectorOfDegree2, Finds0_75InKnotSpanFive) { // NOLINT
  ASSERT_THAT(knot_vector_.GetKnotSpan(ParametricCoordinate{0.75}), Eq(KnotSpan{5}));
}

TEST_F(AKnotVectorOfDegree2, CanCheckIfParametricCoordinateIsEqualLastKnot) { // NOLINT
  ASSERT_THAT(knot_vector_.IsLastKnot(knot_vector_.GetLastKnot()), Eq(true));
}

TEST_F(AKnotVectorOfDegree2, CanCheckIfCoordinateIsNotEqualLastKnot) { // NOLINT
  ASSERT_THAT(knot_vector_.IsLastKnot(ParametricCoordinate{0.9}), Eq(false));
}

//TEST_F(AKnotVectorOfDegree2, CanBeChangedWithAccessOperator) { // NOLINT
//  knot_vector_[0] = ParametricCoordinate{7.0};
//  ASSERT_THAT(knot_vector_[0].Get(), DoubleEq(7.0));
//}

TEST_F(AKnotVectorOfDegree2, FindsParametricCoordinateInKnotVectorRange) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInRange(ParametricCoordinate{0.4}), Eq(true));
}

TEST_F(AKnotVectorOfDegree2, FindsSmallestKnotInKnotVectorRange) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInRange(ParametricCoordinate{0.0}), Eq(true));
}

TEST_F(AKnotVectorOfDegree2, FindsLargestKnotInKnotVectorRange) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInRange(ParametricCoordinate{1.0}), Eq(true));
}

TEST_F(AKnotVectorOfDegree2, DoesNotFindSmallParametricCoordinateInKnotVectorRange) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInRange(ParametricCoordinate{-0.4}), Eq(false));
}

TEST_F(AKnotVectorOfDegree2, DoesNotFindLargeParametricCoordinateInKnotVectorRange) { // NOLINT
  ASSERT_THAT(knot_vector_.IsInRange(ParametricCoordinate{1.5}), Eq(false));
}

TEST_F(AKnotVectorOfDegree2, CanInsertKnot) {  // NOLINT
  baf::KnotVector knot_vector_copy = baf::KnotVector(knot_vector_.begin(), knot_vector_.end());
  knot_vector_copy.InsertKnot(ParametricCoordinate{0.5});
  ASSERT_THAT(knot_vector_copy.GetNumberOfKnots(), knot_vector_.GetNumberOfKnots() + 1);
  ASSERT_THAT(knot_vector_copy.GetKnotSpan(ParametricCoordinate{0.5}).Get(),
              knot_vector_.GetKnotSpan(ParametricCoordinate{0.5}).Get() + 1);
}

TEST_F(AKnotVectorOfDegree2, CanRemoveKnot) {  // NOLINT
  baf::KnotVector knot_vector_copy = baf::KnotVector(knot_vector_.begin(), knot_vector_.end());
  knot_vector_copy.RemoveKnot(ParametricCoordinate{0.5});
  ASSERT_THAT(knot_vector_copy.GetNumberOfKnots(), knot_vector_.GetNumberOfKnots() - 1);
  ASSERT_THAT(knot_vector_copy.GetKnotSpan(ParametricCoordinate{0.5}).Get(),
              knot_vector_.GetKnotSpan(ParametricCoordinate{0.5}).Get() - 1);
}

// TODO(Corinna): commented out the two tests below as the functionality they test is no longer part of class KnotVector
//  need to test this functionality elsewhere --> surface generator.
//TEST_F(AKnotVectorOfDegree2, CanBeAveraged) { // NOLINT
//  std::vector<ParametricCoordinate>
//      coords = {ParametricCoordinate(0.0), ParametricCoordinate(5.0 / 17.0), //ParametricCoordinate(9.0 / 17.0),
//                ParametricCoordinate(14.0 / 17.0), ParametricCoordinate(1.0)};
//  Degree degreeTest(3);
//  int nbControlPoints = 5;
//  baf::KnotVector knot_vector = baf::KnotVector(coords, degreeTest, nbControlPoints);
//  ASSERT_THAT(knot_vector.GetKnot(4).Get(), DoubleEq(28.0 / 51.0));
//  ASSERT_THAT(knot_vector.GetKnot(5).Get(), DoubleEq(1.0));
//  ASSERT_THAT(knot_vector.GetKnot(1).Get(), DoubleEq(0.0));
//}
//
//TEST_F(AKnotVectorOfDegree2, CanBeAveragedLen2) { // NOLINT
//  std::vector<ParametricCoordinate>
//      coords = {ParametricCoordinate(0.0), ParametricCoordinate(10.0 / 17.0), ParametricCoordinate(18.0 / 17.0),
//                ParametricCoordinate(28.0 / 17.0), ParametricCoordinate(2.0)};
//  Degree degreeTest(3);
//  int nbControlPoints = 5;
//  baf::KnotVector knot_vector = baf::KnotVector(coords, degreeTest, nbControlPoints);
//  ASSERT_THAT(knot_vector.GetKnot(4).Get(), DoubleEq(56.0 / 51.0));
//  ASSERT_THAT(knot_vector.GetKnot(5).Get(), DoubleEq(2.0));
//  ASSERT_THAT(knot_vector.GetKnot(1).Get(), DoubleEq(0.0));
//}
