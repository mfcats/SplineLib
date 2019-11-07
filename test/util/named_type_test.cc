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
  TestNamedIntegerType test_named_integer_type_a_{10};
  TestNamedIntegerType test_named_integer_type_b_{-2};
};

TEST_F(TwoNamedIntegerTypes, CanBeCopiedAndComparedIfEqual) {  // NOLINT
  TestNamedIntegerType copy_of_a(test_named_integer_type_a_);
  ASSERT_THAT(copy_of_a, test_named_integer_type_a_);
  TestNamedIntegerType copy_of_b(test_named_integer_type_b_);
  ASSERT_THAT(copy_of_b, test_named_integer_type_b_);
}

TEST_F(TwoNamedIntegerTypes, CanBeAssignedAndComparedIfEqual) {  // NOLINT
  TestNamedIntegerType test_named_integer_type_a_to_be_assigned;
  test_named_integer_type_a_to_be_assigned = test_named_integer_type_a_;
  ASSERT_THAT(test_named_integer_type_a_to_be_assigned, test_named_integer_type_a_);
  TestNamedIntegerType test_named_integer_type_b_to_be_assigned;
  test_named_integer_type_b_to_be_assigned = test_named_integer_type_b_;
  ASSERT_THAT(test_named_integer_type_b_to_be_assigned, test_named_integer_type_b_);
}

TEST_F(TwoNamedIntegerTypes, CanBeIncrementedTo11AndMinus1WithPostIncrementOperator) {  // NOLINT
  ASSERT_THAT((test_named_integer_type_a_++).Get(), Eq(10));
  ASSERT_THAT(test_named_integer_type_a_.Get(), Eq(11));
  ASSERT_THAT((test_named_integer_type_b_++).Get(), Eq(-2));
  ASSERT_THAT(test_named_integer_type_b_.Get(), Eq(-1));
}

TEST_F(TwoNamedIntegerTypes, CanBeIncrementedTo11AndMinus1WithPreIncrementOperator) {  // NOLINT
  ASSERT_THAT((++test_named_integer_type_a_).Get(), Eq(11));
  ASSERT_THAT((++test_named_integer_type_b_).Get(), Eq(-1));
}

TEST_F(TwoNamedIntegerTypes, CanBeDecrementedTo9AndMinus3WithPostDecrementOperator) {  // NOLINT
  ASSERT_THAT((test_named_integer_type_a_--).Get(), Eq(10));
  ASSERT_THAT(test_named_integer_type_a_.Get(), Eq(9));
  ASSERT_THAT((test_named_integer_type_b_--).Get(), Eq(-2));
  ASSERT_THAT(test_named_integer_type_b_.Get(), Eq(-3));
}

TEST_F(TwoNamedIntegerTypes, CanBeDecrementedTo9AndMinus3WithPreDecrementOperator) {  // NOLINT
  ASSERT_THAT((--test_named_integer_type_a_).Get(), Eq(9));
  ASSERT_THAT((--test_named_integer_type_b_).Get(), Eq(-3));
}

TEST_F(TwoNamedIntegerTypes, CanBeModifiedWithNonConstGetMethod) {  // NOLINT
  test_named_integer_type_a_.Get() = -7;
  ASSERT_THAT(test_named_integer_type_a_.Get(), Eq(-7));
  test_named_integer_type_b_.Get() = 0;
  ASSERT_THAT(test_named_integer_type_b_.Get(), Eq(0));
}

TEST_F(TwoNamedIntegerTypes, CanBeComparedIfFirstIsLessThanSecond) {  // NOLINT
  ASSERT_THAT(test_named_integer_type_a_ < test_named_integer_type_b_, false);
  ASSERT_THAT(test_named_integer_type_b_ < test_named_integer_type_a_, true);
}

TEST_F(TwoNamedIntegerTypes, CanBeComparedIfFirstIsLessOrEqualSecond) {  // NOLINT
  ASSERT_THAT(test_named_integer_type_a_ <= test_named_integer_type_b_, false);
  ASSERT_THAT(test_named_integer_type_a_ <= test_named_integer_type_a_, true);
  ASSERT_THAT(test_named_integer_type_b_ <= test_named_integer_type_a_, true);
}

TEST_F(TwoNamedIntegerTypes, CanBeComparedIfFirstIsGreaterThanSecond) {  // NOLINT
  ASSERT_THAT(test_named_integer_type_a_ > test_named_integer_type_b_, true);
  ASSERT_THAT(test_named_integer_type_b_ > test_named_integer_type_a_, false);
}

TEST_F(TwoNamedIntegerTypes, CanBeComparedIfFirstIsGreaterOrEqualSecond) {  // NOLINT
  ASSERT_THAT(test_named_integer_type_a_ >= test_named_integer_type_b_, true);
  ASSERT_THAT(test_named_integer_type_a_ >= test_named_integer_type_a_, true);
  ASSERT_THAT(test_named_integer_type_b_ >= test_named_integer_type_a_, false);
}

TEST_F(TwoNamedIntegerTypes, CanBeAddedTo8) {  // NOLINT
  ASSERT_THAT((test_named_integer_type_a_ + test_named_integer_type_b_).Get(), Eq(8));
  ASSERT_THAT((test_named_integer_type_b_ + test_named_integer_type_a_).Get(), Eq(8));
}

TEST_F(TwoNamedIntegerTypes, CanBeSubtractedToPlusOrMinus12) {  // NOLINT
  ASSERT_THAT((test_named_integer_type_a_ - test_named_integer_type_b_).Get(), Eq(12));
  ASSERT_THAT((test_named_integer_type_b_ - test_named_integer_type_a_).Get(), Eq(-12));
}

class TwoNamedDoubleTypes : public Test {
 public:
  TwoNamedDoubleTypes() = default;

 protected:
  TestNamedDoubleType test_named_double_type_a_{10.01};
  TestNamedDoubleType test_named_double_type_b_{-2.5701};
};

TEST_F(TwoNamedDoubleTypes, CanBeCopiedAndComparedIfEqual) {  // NOLINT
  TestNamedDoubleType copy_of_a(test_named_double_type_a_);
  ASSERT_THAT(copy_of_a, test_named_double_type_a_);
  TestNamedDoubleType copy_of_b(test_named_double_type_b_);
  ASSERT_THAT(copy_of_b, test_named_double_type_b_);
}

TEST_F(TwoNamedDoubleTypes, CanBeAssignedAndComparedIfEqual) {  // NOLINT
  TestNamedDoubleType test_named_double_type_a_to_be_assigned;
  test_named_double_type_a_to_be_assigned = test_named_double_type_a_;
  ASSERT_THAT(test_named_double_type_a_to_be_assigned, test_named_double_type_a_);
  TestNamedDoubleType test_named_double_type_b_to_be_assigned;
  test_named_double_type_b_to_be_assigned = test_named_double_type_b_;
  ASSERT_THAT(test_named_double_type_b_to_be_assigned, test_named_double_type_b_);
}

TEST_F(TwoNamedDoubleTypes, CanBeIncrementedTo11_01AndMinus1_5701WithPostIncrementOperator) {  // NOLINT
  ASSERT_THAT((test_named_double_type_a_++).Get(), DoubleEq(10.01));
  ASSERT_THAT(test_named_double_type_a_.Get(), DoubleEq(11.01));
  ASSERT_THAT((test_named_double_type_b_++).Get(), DoubleEq(-2.5701));
  ASSERT_THAT(test_named_double_type_b_.Get(), DoubleEq(-1.5701));
}

TEST_F(TwoNamedDoubleTypes, CanBeIncrementedTo11_01AndMinus1_5701WithPreIncrementOperator) {  // NOLINT
  ASSERT_THAT((++test_named_double_type_a_).Get(), DoubleEq(11.01));
  ASSERT_THAT((++test_named_double_type_b_).Get(), DoubleEq(-1.5701));
}

TEST_F(TwoNamedDoubleTypes, CanBeDecrementedTo9_01AndMinus3_5701WithPostDecrementOperator) {  // NOLINT
  ASSERT_THAT((test_named_double_type_a_--).Get(), DoubleEq(10.01));
  ASSERT_THAT(test_named_double_type_a_.Get(), DoubleEq(9.01));
  ASSERT_THAT((test_named_double_type_b_--).Get(), DoubleEq(-2.5701));
  ASSERT_THAT(test_named_double_type_b_.Get(), DoubleEq(-3.5701));
}

TEST_F(TwoNamedDoubleTypes, CanBeDecrementedTo9_01AndMinus3_5701WithPreDecrementOperator) {  // NOLINT
  ASSERT_THAT((--test_named_double_type_a_).Get(), DoubleEq(9.01));
  ASSERT_THAT((--test_named_double_type_b_).Get(), DoubleEq(-3.5701));
}

TEST_F(TwoNamedDoubleTypes, CanBeModifiedWithNonConstGetMethod) {  // NOLINT
  test_named_double_type_a_.Get() = -7.909;
  ASSERT_THAT(test_named_double_type_a_.Get(), DoubleEq(-7.909));
  test_named_double_type_b_.Get() = 0.0;
  ASSERT_THAT(test_named_double_type_b_.Get(), DoubleEq(0.0));
}

TEST_F(TwoNamedDoubleTypes, CanBeComparedIfFirstIsLessThanSecond) {  // NOLINT
  ASSERT_THAT(test_named_double_type_a_ < test_named_double_type_b_, false);
  ASSERT_THAT(test_named_double_type_b_ < test_named_double_type_a_, true);
}

TEST_F(TwoNamedDoubleTypes, CanBeComparedIfFirstIsLessOrEqualSecond) {  // NOLINT
  ASSERT_THAT(test_named_double_type_a_ <= test_named_double_type_b_, false);
  ASSERT_THAT(test_named_double_type_a_ <= test_named_double_type_a_, true);
  ASSERT_THAT(test_named_double_type_b_ <= test_named_double_type_a_, true);
}

TEST_F(TwoNamedDoubleTypes, CanBeComparedIfFirstIsGreaterThanSecond) {  // NOLINT
  ASSERT_THAT(test_named_double_type_a_ > test_named_double_type_b_, true);
  ASSERT_THAT(test_named_double_type_b_ > test_named_double_type_a_, false);
}

TEST_F(TwoNamedDoubleTypes, CanBeComparedIfFirstIsGreaterOrEqualSecond) {  // NOLINT
  ASSERT_THAT(test_named_double_type_a_ >= test_named_double_type_b_, true);
  ASSERT_THAT(test_named_double_type_a_ >= test_named_double_type_a_, true);
  ASSERT_THAT(test_named_double_type_b_ >= test_named_double_type_a_, false);
}

TEST_F(TwoNamedDoubleTypes, CanBeAddedTo7_4399) {  // NOLINT
  ASSERT_THAT((test_named_double_type_a_ + test_named_double_type_b_).Get(), DoubleEq(7.4399));
  ASSERT_THAT((test_named_double_type_b_ + test_named_double_type_a_).Get(), DoubleEq(7.4399));
}

TEST_F(TwoNamedDoubleTypes, CanBeSubtractedToPlusOrMinus12_5801) {  // NOLINT
  ASSERT_THAT((test_named_double_type_a_ - test_named_double_type_b_).Get(), DoubleEq(12.5801));
  ASSERT_THAT((test_named_double_type_b_ - test_named_double_type_a_).Get(), DoubleEq(-12.5801));
}
