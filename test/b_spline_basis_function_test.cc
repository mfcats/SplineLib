/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University
This file is part of SplineLib.
SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.
SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <memory>
#include <cstdint>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "b_spline_basis_function.h"
#include "named_type.h"

using ParamCoord = util::NamedType<double, struct ParamCoordParameter>;

using std::make_shared;
using std::vector;

using testing::DoubleEq;
using testing::Test;
using ::testing::AtLeast;
using ::testing::Return;

class MockKnotVector : public baf::KnotVector {
  public:
   MOCK_CONST_METHOD0(GetNumberOfKnots, size_t());
   MOCK_CONST_METHOD1(IsInKnotVectorRange, bool(ParamCoord));
   MOCK_CONST_METHOD1(GetKnotSpan, int64_t(ParamCoord));
   MOCK_CONST_METHOD1(GetKnot, ParamCoord(size_t));
 };

// Test basis function N_{0,1} from NURBS book example 2.1
class BasisFunctionEx21N01 : public Test {
 public:
  BasisFunctionEx21N01() :
      knot_vector_ptr(std::make_shared<MockKnotVector>()),
      basis_function_(knot_vector_ptr, 1, 0) {}

 protected:
  std::shared_ptr<MockKnotVector> knot_vector_ptr;
  baf::BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N01, Mocking) {
  EXPECT_CALL(*knot_vector_ptr, GetNumberOfKnots())
      .WillOnce(Return(3));
  ASSERT_EQ(basis_function_.testMock(), 3);
}

TEST_F(BasisFunctionEx21N01, IsZeroAt0_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.0}))
      .Times(1)
      .WillOnce(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.0}))
      .Times(1)
      .WillOnce(Return(int64_t(2)));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N01, IsZeroAt0_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.5}))
      .Times(1)
      .WillOnce(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.5}))
      .Times(1)
      .WillOnce(Return(2));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N01, IsZeroAt1_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.0}))
      .Times(1)
      .WillOnce(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{1.0}))
      .Times(1)
      .WillOnce(Return(2));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N01, IsZeroAt1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.5}))
      .Times(1)
      .WillOnce(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N01, IsZeroAtMinus1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{-1.5}))
      .Times(1)
      .WillOnce(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{1,1} from NURBS book example 2.1
class BasisFunctionEx21N11 : public Test {
 public:
  BasisFunctionEx21N11() :
      knot_vector_ptr(std::make_shared<MockKnotVector>()),
      basis_function_(knot_vector_ptr, 1, 1) {}

 protected:
  std::shared_ptr<MockKnotVector> knot_vector_ptr;
  baf::BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N11, IsOneAt0_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.0}))
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.0}))
      .Times(3)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(1))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.0}), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx21N11, IsHalfAt0_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.5}))
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.5}))
      .Times(3)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(1))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.5}), DoubleEq(0.5));
}

TEST_F(BasisFunctionEx21N11, IsZeroAt1_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.0}))
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{1.0}))
      .Times(3)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(1))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N11, IsZeroAt1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.5}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N11, IsZeroAtMinus1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{-1.5}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{-1.5}), DoubleEq(0.0));
}


// Test basis function N_{2,1} from NURBS book example 2.1
class BasisFunctionEx21N21 : public Test {
 public:
  BasisFunctionEx21N21() :
      knot_vector_ptr(std::make_shared<MockKnotVector>()),
      basis_function_(knot_vector_ptr, 1, 2) {}

 protected:
  std::shared_ptr<MockKnotVector> knot_vector_ptr;
  baf::BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N21, IsZeroAt0_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.0}))
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.0}))
      .Times(3)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N21, IsHalfAt0_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.5}))
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.5}))
      .Times(3)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.5}), DoubleEq(0.5));
}

TEST_F(BasisFunctionEx21N21, IsOneAt1_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.0}))
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{1.0}))
      .Times(3)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.0}), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx21N21, IsZeroAt1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.5}))
      .Times(1)
      .WillOnce(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N21, IsZeroAtMinus1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{-1.5}))
      .Times(1)
      .WillOnce(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{3,1} from NURBS book example 2.1
class BasisFunctionEx21N31 : public Test {
 public:
  BasisFunctionEx21N31() :
      knot_vector_ptr(std::make_shared<MockKnotVector>()),
      basis_function_(knot_vector_ptr, 1, 3) {}

 protected:
  std::shared_ptr<MockKnotVector> knot_vector_ptr;
  baf::BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N31, IsZeroAt0_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.0}))
      .Times(1)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.0}))
      .Times(1)
      .WillRepeatedly(Return(int64_t(2)));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N31, IsZeroAt0_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.5}))
      .Times(1)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.5}))
      .Times(1)
      .WillRepeatedly(Return(int64_t(2)));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N31, IsZeroAt1_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.0}))
      .Times(1)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{1.0}))
      .Times(1)
      .WillRepeatedly(Return(int64_t(2)));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N31, IsZeroAt1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.5}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N31, IsZeroAtMinus1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{-1.5}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{0,2} from NURBS book example 2.1
class BasisFunctionEx21N02 : public Test {
 public:
  BasisFunctionEx21N02() :
      knot_vector_ptr(std::make_shared<MockKnotVector>()),
      basis_function_(knot_vector_ptr, 2, 0) {}

 protected:
  std::shared_ptr<MockKnotVector> knot_vector_ptr;
  baf::BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N02, IsZeroAt0_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.0}))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.0}))
      .Times(5)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(0))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(1))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(5)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.0}), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx21N02, Is0_25At0_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.5}))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.5}))
      .Times(5)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(0))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(1))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(5)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.5}), DoubleEq(0.25));
}

TEST_F(BasisFunctionEx21N02, IsZeroAt1_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.0}))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{1.0}))
      .Times(5)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(0))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(1))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(5)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N02, IsZeroAt1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.5}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N02, IsZeroAtMinus1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{-1.5}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{1,2} from NURBS book example 2.1
class BasisFunctionEx21N12 : public Test {
 public:
  BasisFunctionEx21N12() :
      knot_vector_ptr(std::make_shared<MockKnotVector>()),
      basis_function_(knot_vector_ptr, 2, 1) {}

 protected:
  std::shared_ptr<MockKnotVector> knot_vector_ptr;
  baf::BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N12, IsZeroAt0_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.0}))
      .Times(7)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.0}))
      .Times(7)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(1))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(7)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(7)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N12, Is0_5At0_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.5}))
      .Times(7)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.5}))
      .Times(7)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(1))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(7)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(7)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.5}), DoubleEq(0.5));
}

TEST_F(BasisFunctionEx21N12, IsZeroAt1_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.0}))
      .Times(7)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{1.0}))
      .Times(7)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(1))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(7)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(7)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N12, IsZeroAt1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.5}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N12, IsZeroAtMinus1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{-1.5}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{2,2} from NURBS book example 2.1
class BasisFunctionEx21N22 : public Test {
 public:
  BasisFunctionEx21N22() :
      knot_vector_ptr(std::make_shared<MockKnotVector>()),
      basis_function_(knot_vector_ptr, 2, 2) {}

 protected:
  std::shared_ptr<MockKnotVector> knot_vector_ptr;
  baf::BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N22, IsZeroAt0_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.0}))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.0}))
      .Times(5)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(5)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(5))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N22, Is0_25At0_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.5}))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.5}))
      .Times(5)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(5)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(5))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.5}), DoubleEq(0.25));
}

TEST_F(BasisFunctionEx21N22, IsOneAt1_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.0}))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{1.0}))
      .Times(5)
      .WillRepeatedly(Return(int64_t(2)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(2))
      .Times(5)
      .WillRepeatedly(Return(ParamCoord{0.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(5))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{1.0}));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.0}), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx21N22, IsZeroAt1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.5}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N22, IsZeroAtMinus1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{-1.5}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{0,1} from NURBS book example 2.2
class BasisFunctionEx22N01 : public Test {
 public:
  BasisFunctionEx22N01() : knot_vector_ptr(std::make_shared<MockKnotVector>()),
                           basis_function_(knot_vector_ptr, 1, 0) {}

 protected:
  std::shared_ptr<MockKnotVector> knot_vector_ptr;
  baf::BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt0_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.0}))
      .Times(2)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.0}))
      .Times(2)
      .WillRepeatedly(Return(int64_t(2)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, 0),
              DoubleEq(basis_function_.Evaluate(ParamCoord{0.0})));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.5}))
      .Times(2)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{1.5}))
      .Times(2)
      .WillRepeatedly(Return(int64_t(3)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1.5}, 0),
              DoubleEq(basis_function_.Evaluate(ParamCoord{1.5})));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt2_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{2.0}))
      .Times(2)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{2.0}))
      .Times(2)
      .WillRepeatedly(Return(int64_t(3)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{2.0}, 0),
              DoubleEq(basis_function_.Evaluate(ParamCoord{2.0})));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt4_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{4.0}))
      .Times(2)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{4.0}))
      .Times(2)
      .WillRepeatedly(Return(int64_t(3)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{4.0}, 0),
              DoubleEq(basis_function_.Evaluate(ParamCoord{4.0})));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt5_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{5}))
      .Times(2)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{5}))
      .Times(2)
      .WillRepeatedly(Return(int64_t(3)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{5.0}, 0),
              DoubleEq(basis_function_.Evaluate(ParamCoord{5.0})));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt6_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{6}))
      .Times(2)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{6}))
      .Times(2)
      .WillRepeatedly(Return(int64_t(3)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{6.0}, 0),
              DoubleEq(basis_function_.Evaluate(ParamCoord{6.0})));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAtMinus0_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{-0.5}))
      .Times(2)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{-0.5}, 0),
              DoubleEq(basis_function_.Evaluate(ParamCoord{-0.5})));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt0_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.0}))
      .Times(1)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.0}))
      .Times(1)
      .WillRepeatedly(Return(int64_t(2)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0.0}, 1), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.5}))
      .Times(1)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{1.5}))
      .Times(1)
      .WillRepeatedly(Return(int64_t(3)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1.5}, 1), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt2_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{2.0}))
      .Times(1)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{2.0}))
      .Times(1)
      .WillRepeatedly(Return(int64_t(4)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{2.0}, 1), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt4_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{4.0}))
      .Times(1)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{4.0}))
      .Times(1)
      .WillRepeatedly(Return(int64_t(7)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{4.0}, 1), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt5_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{5}))
      .Times(1)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{5}))
      .Times(1)
      .WillRepeatedly(Return(int64_t(9)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{5.0}, 1), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt6_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{6}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{6.0}, 1), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAtMinus0_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{-0.5}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{-0.5}, 1), DoubleEq(0.0));
}

// Test basis function derivative N_{3,1} from NURBS book example 2.2
class BasisFunctionEx22N13 : public Test {
 public:
  BasisFunctionEx22N13() : knot_vector_ptr(std::make_shared<MockKnotVector>()),
                           basis_function_(knot_vector_ptr, 1, 3) {}

 protected:
  std::shared_ptr<MockKnotVector> knot_vector_ptr;
  baf::BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqualZeroAt0_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.0}))
       .Times(1)
       .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.0}))
       .Times(1)
       .WillRepeatedly(Return(int64_t(2)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0.0}, 1), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqualZeroAt0_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{0.5}))
      .Times(1)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{0.5}))
      .Times(1)
      .WillRepeatedly(Return(int64_t(2)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0.5}, 1), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqual1At1_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.0}))
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{1.0}))
      .Times(3)
      .WillRepeatedly(Return(int64_t(3)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{1}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{2.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(5))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{3.0}));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1.0}, 1), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqual1At1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.5}))
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{1.5}))
      .Times(3)
      .WillRepeatedly(Return(int64_t(3)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{1}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{2.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(5))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{3.0}));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1.5}, 1), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqualMinus1At2_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{2}))
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{2}))
      .Times(3)
      .WillRepeatedly(Return(int64_t(4)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{1}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{2.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(5))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{3.0}));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{2.0}, 1), DoubleEq(-1.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqualMinus1At2_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{2.5}))
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{2.5}))
      .Times(3)
      .WillRepeatedly(Return(int64_t(4)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{1}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{2.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(5))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{3.0}));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{2.5}, 1), DoubleEq(-1.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqual0At3_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{3}))
      .Times(1)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{3}))
      .Times(1)
      .WillRepeatedly(Return(int64_t(5)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{3.0}, 1), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N13, FourthDerevitveIsEqual0At1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.5}))
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{1.5}))
      .Times(3)
      .WillRepeatedly(Return(int64_t(3)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(3))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{1}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(4))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{2.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(5))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{3.0}));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1.5}, 4), DoubleEq(0.0));
}

// Test basis function derivative N_{6,1} from NURBS book example 2.2
class BasisFunctionEx22N61 : public Test {
 public:
  BasisFunctionEx22N61() : knot_vector_ptr(std::make_shared<MockKnotVector>()),
                           basis_function_(knot_vector_ptr, 1, 6) {}

 protected:
  std::shared_ptr<MockKnotVector> knot_vector_ptr;
  baf::BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx22N61, FirstDerevitveIsEqualMinus1At4_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{4.0}))
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{4.0}))
      .Times(3)
      .WillRepeatedly(Return(int64_t(7)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(6))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{1}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(7))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{2.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(8))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{3.0}));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{4.0}, 1), DoubleEq(-1.0));
}

TEST_F(BasisFunctionEx22N61, SecondDerevitveIsEqual0At4_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{4.5}))
      .Times(3)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{4.5}))
      .Times(3)
      .WillRepeatedly(Return(int64_t(7)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(6))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{1}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(7))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{2.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(8))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{3.0}));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{4.5}, 2), DoubleEq(0.0));
}

// Test basis function derivative N_{7,2} from NURBS book example 2.2
class BasisFunctionEx22N72 : public Test {
 public:
  BasisFunctionEx22N72() : knot_vector_ptr(std::make_shared<MockKnotVector>()),
                           basis_function_(knot_vector_ptr, 2, 7) {}

 protected:
  std::shared_ptr<MockKnotVector> knot_vector_ptr;
  baf::BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual0At4_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{4}))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{4}))
      .Times(5)
      .WillRepeatedly(Return(int64_t(7)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(7))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{4}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(8))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(9))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(10))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{4.0}, 1), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, SecondDerevitveIsEqual2At4_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{4}))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{4}))
      .Times(5)
      .WillRepeatedly(Return(int64_t(7)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(7))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{4}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(8))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(9))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(10))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{4.0}, 2), DoubleEq(2.0));
}

TEST_F(BasisFunctionEx22N72, ThirdDerevitveIsEqual0At4_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{4}))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{4}))
      .Times(5)
      .WillRepeatedly(Return(int64_t(7)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(7))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{4}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(8))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(9))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(10))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{4.0}, 3), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual1At4_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{4.5}))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{4.5}))
      .Times(5)
      .WillRepeatedly(Return(int64_t(7)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(7))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{4}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(8))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(9))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(10))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{4.5}, 1), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx22N72, SecondDerevitveIsEqual2At4_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{4.5}))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{4.5}))
      .Times(5)
      .WillRepeatedly(Return(int64_t(7)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(7))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{4}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(8))
      .Times(4)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(9))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(10))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{4.5}, 2), DoubleEq(2.0));
}

TEST_F(BasisFunctionEx22N72, ThirdDerevitveIsEqual0At4_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{4.5}))
      .Times(5)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{4.5}))
      .Times(5)
      .WillRepeatedly(Return(int64_t(7)));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(7))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{4}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(8))
      .Times(3)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(9))
      .Times(2)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  EXPECT_CALL(*knot_vector_ptr, GetKnot(10))
      .Times(1)
      .WillRepeatedly(Return(ParamCoord{5.0}));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{4.5}, 3), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual2At5_0) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{5}))
      .Times(1)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{5}))
      .Times(1)
      .WillRepeatedly(Return(int64_t(10)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{5.0}, 1), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual0At5_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{5.5}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{5.5}, 1), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual0At1_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{1.5}))
      .Times(1)
      .WillRepeatedly(Return(true));
  EXPECT_CALL(*knot_vector_ptr, GetKnotSpan(ParamCoord{1.5}))
      .Times(1)
      .WillRepeatedly(Return(int64_t(3)));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1.5}, 1), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual0AtMinus5_5) {
  EXPECT_CALL(*knot_vector_ptr, IsInKnotVectorRange(ParamCoord{-5.5}))
      .Times(1)
      .WillRepeatedly(Return(false));
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{-5.5}, 1), DoubleEq(0.0));
}
