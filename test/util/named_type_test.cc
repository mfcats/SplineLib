/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/util/named_type.h"

using testing::DoubleEq;
using testing::Eq;
using testing::Test;

using namespace splinelib::src;
using TestNamedIntegerType  = util::NamedType<int, struct TestNamedIntegerTypeName>;
using TestNamedDoubleType  = util::NamedType<double, struct TestNamedDoubleTypeName>;

class TwoNamedIntegerTypes : public Test {
 public:
  TwoNamedIntegerTypes() = default;

 protected:
  TestNamedIntegerType named_integer_type_a_{10};
  TestNamedIntegerType named_integer_type_b_{-2};
};

TEST_F(TwoNamedIntegerTypes, CanBeCopiedAndComparedIfEqual) {  // NOLINT
  TestNamedIntegerType copied_named_integer_type_a = named_integer_type_a_;
  ASSERT_THAT(copied_named_integer_type_a, named_integer_type_a_);
}

TEST_F(TwoNamedIntegerTypes, CanBeCreatedWithMoveConstructor) {  // NOLINT
  TestNamedIntegerType moved_named_integer_type_a(std::move(named_integer_type_a_));
  ASSERT_THAT(moved_named_integer_type_a.Get(), Eq(10));
}

TEST_F(TwoNamedIntegerTypes, CanBeAssignedAndComparedIfEqual) {  // NOLINT
  TestNamedIntegerType named_integer_type_a_to_assign{};
  named_integer_type_a_to_assign = named_integer_type_a_;
  ASSERT_THAT(named_integer_type_a_to_assign, named_integer_type_a_);
}

TEST_F(TwoNamedIntegerTypes, CanBeMoveAssignedAndComparedIfEqual) {  // NOLINT
  TestNamedIntegerType named_integer_type_a_to_move_assign{};
  named_integer_type_a_to_move_assign = std::move(named_integer_type_a_);
  ASSERT_THAT(named_integer_type_a_to_move_assign.Get(), Eq(10));
}

TEST_F(TwoNamedIntegerTypes, CanBeIncrementedTo11AndMinus1WithPostIncrementOperator) {  // NOLINT
  ASSERT_THAT((named_integer_type_a_++).Get(), Eq(10));
  ASSERT_THAT(named_integer_type_a_.Get(), Eq(11));
  ASSERT_THAT((named_integer_type_b_++).Get(), Eq(-2));
  ASSERT_THAT(named_integer_type_b_.Get(), Eq(-1));
}

TEST_F(TwoNamedIntegerTypes, CanBeIncrementedTo11AndMinus1WithPreIncrementOperator) {  // NOLINT
  ASSERT_THAT((++named_integer_type_a_).Get(), Eq(11));
  ASSERT_THAT((++named_integer_type_b_).Get(), Eq(-1));
}

TEST_F(TwoNamedIntegerTypes, CanBeDecrementedTo9AndMinus3WithPostDecrementOperator) {  // NOLINT
  ASSERT_THAT((named_integer_type_a_--).Get(), Eq(10));
  ASSERT_THAT(named_integer_type_a_.Get(), Eq(9));
  ASSERT_THAT((named_integer_type_b_--).Get(), Eq(-2));
  ASSERT_THAT(named_integer_type_b_.Get(), Eq(-3));
}

TEST_F(TwoNamedIntegerTypes, CanBeDecrementedTo9AndMinus3WithPreDecrementOperator) {  // NOLINT
  ASSERT_THAT((--named_integer_type_a_).Get(), Eq(9));
  ASSERT_THAT((--named_integer_type_b_).Get(), Eq(-3));
}

TEST_F(TwoNamedIntegerTypes, CanBeModifiedWithNonConstGetMethod) {  // NOLINT
  named_integer_type_a_.Get() = -7;
  ASSERT_THAT(named_integer_type_a_.Get(), Eq(-7));
  named_integer_type_b_.Get() = 0;
  ASSERT_THAT(named_integer_type_b_.Get(), Eq(0));
}

TEST_F(TwoNamedIntegerTypes, CanBeComparedIfFirstIsLessThanSecond) {  // NOLINT
  ASSERT_THAT(named_integer_type_a_ < named_integer_type_b_, false);
  ASSERT_THAT(named_integer_type_b_ < named_integer_type_a_, true);
}

TEST_F(TwoNamedIntegerTypes, CanBeComparedIfFirstIsLessOrEqualSecond) {  // NOLINT
  ASSERT_THAT(named_integer_type_a_ <= named_integer_type_b_, false);
  ASSERT_THAT(named_integer_type_a_ <= named_integer_type_a_, true);
  ASSERT_THAT(named_integer_type_b_ <= named_integer_type_a_, true);
}

TEST_F(TwoNamedIntegerTypes, CanBeComparedIfFirstIsGreaterThanSecond) {  // NOLINT
  ASSERT_THAT(named_integer_type_a_ > named_integer_type_b_, true);
  ASSERT_THAT(named_integer_type_b_ > named_integer_type_a_, false);
}

TEST_F(TwoNamedIntegerTypes, CanBeComparedIfFirstIsGreaterOrEqualSecond) {  // NOLINT
  ASSERT_THAT(named_integer_type_a_ >= named_integer_type_b_, true);
  ASSERT_THAT(named_integer_type_a_ >= named_integer_type_a_, true);
  ASSERT_THAT(named_integer_type_b_ >= named_integer_type_a_, false);
}

TEST_F(TwoNamedIntegerTypes, CanBeAddedTo8) {  // NOLINT
  ASSERT_THAT((named_integer_type_a_ + named_integer_type_b_).Get(), Eq(8));
  ASSERT_THAT((named_integer_type_b_ + named_integer_type_a_).Get(), Eq(8));
}

TEST_F(TwoNamedIntegerTypes, CanBeSubtractedToPlusOrMinus12) {  // NOLINT
  ASSERT_THAT((named_integer_type_a_ - named_integer_type_b_).Get(), Eq(12));
  ASSERT_THAT((named_integer_type_b_ - named_integer_type_a_).Get(), Eq(-12));
}

TEST_F(TwoNamedIntegerTypes, CanBeMulitpliedByPlusAndMinus2To20And4) {  // NOLINT
  ASSERT_THAT((2 * named_integer_type_a_).Get(), Eq(20));
  ASSERT_THAT((named_integer_type_a_ * 2).Get(), Eq(20));
  ASSERT_THAT((-2 * named_integer_type_b_).Get(), Eq(4));
  ASSERT_THAT((named_integer_type_b_ * (-2)).Get(), Eq(4));
}

class TwoNamedDoubleTypes : public Test {
 public:
  TwoNamedDoubleTypes() = default;

 protected:
  TestNamedDoubleType named_double_type_a_{10.01};
  TestNamedDoubleType named_double_type_b_{-2.5701};
};

TEST_F(TwoNamedDoubleTypes, CanBeCopiedAndComparedIfEqual) {  // NOLINT
  TestNamedDoubleType copied_named_double_type_a = named_double_type_a_;
  ASSERT_THAT(copied_named_double_type_a, named_double_type_a_);
}

TEST_F(TwoNamedDoubleTypes, CanBeCreatedWithMoveConstructor) {  // NOLINT
  TestNamedDoubleType moved_named_double_type_a(std::move(named_double_type_a_));
  ASSERT_THAT(moved_named_double_type_a.Get(), DoubleEq(10.01));
}

TEST_F(TwoNamedDoubleTypes, CanBeAssignedAndComparedIfEqual) {  // NOLINT
  TestNamedDoubleType named_double_type_a_to_assign{};
  named_double_type_a_to_assign = named_double_type_a_;
  ASSERT_THAT(named_double_type_a_to_assign, named_double_type_a_);
}

TEST_F(TwoNamedDoubleTypes, CanBeMoveAssignedAndComparedIfEqual) {  // NOLINT
  TestNamedDoubleType named_double_type_a_to_move_assign{};
  named_double_type_a_to_move_assign = std::move(named_double_type_a_);
  ASSERT_THAT(named_double_type_a_to_move_assign.Get(), DoubleEq(10.01));
}

TEST_F(TwoNamedDoubleTypes, CanBeIncrementedTo11_01AndMinus1_5701WithPostIncrementOperator) {  // NOLINT
  ASSERT_THAT((named_double_type_a_++).Get(), DoubleEq(10.01));
  ASSERT_THAT(named_double_type_a_.Get(), DoubleEq(11.01));
  ASSERT_THAT((named_double_type_b_++).Get(), DoubleEq(-2.5701));
  ASSERT_THAT(named_double_type_b_.Get(), DoubleEq(-1.5701));
}

TEST_F(TwoNamedDoubleTypes, CanBeIncrementedTo11_01AndMinus1_5701WithPreIncrementOperator) {  // NOLINT
  ASSERT_THAT((++named_double_type_a_).Get(), DoubleEq(11.01));
  ASSERT_THAT((++named_double_type_b_).Get(), DoubleEq(-1.5701));
}

TEST_F(TwoNamedDoubleTypes, CanBeDecrementedTo9_01AndMinus3_5701WithPostDecrementOperator) {  // NOLINT
  ASSERT_THAT((named_double_type_a_--).Get(), DoubleEq(10.01));
  ASSERT_THAT(named_double_type_a_.Get(), DoubleEq(9.01));
  ASSERT_THAT((named_double_type_b_--).Get(), DoubleEq(-2.5701));
  ASSERT_THAT(named_double_type_b_.Get(), DoubleEq(-3.5701));
}

TEST_F(TwoNamedDoubleTypes, CanBeDecrementedTo9_01AndMinus3_5701WithPreDecrementOperator) {  // NOLINT
  ASSERT_THAT((--named_double_type_a_).Get(), DoubleEq(9.01));
  ASSERT_THAT((--named_double_type_b_).Get(), DoubleEq(-3.5701));
}

TEST_F(TwoNamedDoubleTypes, CanBeModifiedWithNonConstGetMethod) {  // NOLINT
  named_double_type_a_.Get() = -7.909;
  ASSERT_THAT(named_double_type_a_.Get(), DoubleEq(-7.909));
  named_double_type_b_.Get() = 0.0;
  ASSERT_THAT(named_double_type_b_.Get(), DoubleEq(0.0));
}

TEST_F(TwoNamedDoubleTypes, CanBeComparedIfFirstIsLessThanSecond) {  // NOLINT
  ASSERT_THAT(named_double_type_a_ < named_double_type_b_, false);
  ASSERT_THAT(named_double_type_b_ < named_double_type_a_, true);
}

TEST_F(TwoNamedDoubleTypes, CanBeComparedIfFirstIsLessOrEqualSecond) {  // NOLINT
  ASSERT_THAT(named_double_type_a_ <= named_double_type_b_, false);
  ASSERT_THAT(named_double_type_a_ <= named_double_type_a_, true);
  ASSERT_THAT(named_double_type_b_ <= named_double_type_a_, true);
}

TEST_F(TwoNamedDoubleTypes, CanBeComparedIfFirstIsGreaterThanSecond) {  // NOLINT
  ASSERT_THAT(named_double_type_a_ > named_double_type_b_, true);
  ASSERT_THAT(named_double_type_b_ > named_double_type_a_, false);
}

TEST_F(TwoNamedDoubleTypes, CanBeComparedIfFirstIsGreaterOrEqualSecond) {  // NOLINT
  ASSERT_THAT(named_double_type_a_ >= named_double_type_b_, true);
  ASSERT_THAT(named_double_type_a_ >= named_double_type_a_, true);
  ASSERT_THAT(named_double_type_b_ >= named_double_type_a_, false);
}

TEST_F(TwoNamedDoubleTypes, CanBeAddedTo7_4399) {  // NOLINT
  ASSERT_THAT((named_double_type_a_ + named_double_type_b_).Get(), DoubleEq(7.4399));
  ASSERT_THAT((named_double_type_b_ + named_double_type_a_).Get(), DoubleEq(7.4399));
}

TEST_F(TwoNamedDoubleTypes, CanBeSubtractedToPlusOrMinus12_5801) {  // NOLINT
  ASSERT_THAT((named_double_type_a_ - named_double_type_b_).Get(), DoubleEq(12.5801));
  ASSERT_THAT((named_double_type_b_ - named_double_type_a_).Get(), DoubleEq(-12.5801));
}

TEST_F(TwoNamedDoubleTypes, CanBeMulitpliedByPlusAndMinus1_5To15_015And3_85515) {  // NOLINT
  ASSERT_THAT((1.5 * named_double_type_a_).Get(), Eq(15.015));
  ASSERT_THAT((named_double_type_a_ * 1.5).Get(), Eq(15.015));
  ASSERT_THAT((-1.5 * named_double_type_b_).Get(), Eq(3.85515));
  ASSERT_THAT((named_double_type_b_ * (-1.5)).Get(), Eq(3.85515));
}
