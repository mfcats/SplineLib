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

#include "gmock/gmock.h"

#include "b_spline_basis_function.h"

using std::make_shared;
using std::vector;

using testing::DoubleEq;
using testing::Test;

// Test basis function N_{0,1} from NURBS book example 2.1
class BasisFunctionEx21N01 : public Test {
 public:
  BasisFunctionEx21N01() :
      knot_vector_(vector<double>({0, 0, 0, 1, 1, 1})),
      basis_function_(knot_vector_, 1, 0) {}

 protected:
  KnotVector knot_vector_;
  BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N01, IsZeroAt0_0) {
  ASSERT_THAT(basis_function_.Evaluate(0.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N01, IsZeroAt0_5) {
  ASSERT_THAT(basis_function_.Evaluate(0.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N01, IsZeroAt1_0) {
  ASSERT_THAT(basis_function_.Evaluate(1.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N01, IsZeroAt1_5) {
  ASSERT_THAT(basis_function_.Evaluate(1.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N01, IsZeroAtMinus1_5) {
  ASSERT_THAT(basis_function_.Evaluate(-1.5), DoubleEq(0.0));
}

// Test basis function N_{1,1} from NURBS book example 2.1
class BasisFunctionEx21N11 : public Test {
 public:
  BasisFunctionEx21N11() :
      knot_vector_(vector<double>({0, 0, 0, 1, 1, 1})),
      basis_function_(knot_vector_, 1, 1) {}

 protected:
  KnotVector knot_vector_;
  BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N11, IsZeroAt0_0) {
  ASSERT_THAT(basis_function_.Evaluate(0.0), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx21N11, IsZeroAt0_5) {
  ASSERT_THAT(basis_function_.Evaluate(0.5), DoubleEq(0.5));
}

TEST_F(BasisFunctionEx21N11, IsZeroAt1_0) {
  ASSERT_THAT(basis_function_.Evaluate(1.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N11, IsZeroAt1_5) {
  ASSERT_THAT(basis_function_.Evaluate(1.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N11, IsZeroAtMinus1_5) {
  ASSERT_THAT(basis_function_.Evaluate(-1.5), DoubleEq(0.0));
}

// Test basis function N_{2,1} from NURBS book example 2.1
class BasisFunctionEx21N21 : public Test {
 public:
  BasisFunctionEx21N21() :
      knot_vector_(vector<double>({0, 0, 0, 1, 1, 1})),
      basis_function_(knot_vector_, 1, 2) {}

 protected:
  KnotVector knot_vector_;
  BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N21, IsZeroAt0_0) {
  ASSERT_THAT(basis_function_.Evaluate(0.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N21, IsZeroAt0_5) {
  ASSERT_THAT(basis_function_.Evaluate(0.5), DoubleEq(0.5));
}

TEST_F(BasisFunctionEx21N21, IsOneAt1_0) {
  ASSERT_THAT(basis_function_.Evaluate(1.0), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx21N21, IsZeroAt1_5) {
  ASSERT_THAT(basis_function_.Evaluate(1.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N21, IsZeroAtMinus1_5) {
  ASSERT_THAT(basis_function_.Evaluate(-1.5), DoubleEq(0.0));
}

// Test basis function N_{3,1} from NURBS book example 2.1
class BasisFunctionEx21N31 : public Test {
 public:
  BasisFunctionEx21N31() :
      knot_vector_(vector<double>({0, 0, 0, 1, 1, 1})),
      basis_function_(knot_vector_, 1, 3) {}

 protected:
  KnotVector knot_vector_;
  BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N31, IsZeroAt0_0) {
  ASSERT_THAT(basis_function_.Evaluate(0.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N31, IsZeroAt0_5) {
  ASSERT_THAT(basis_function_.Evaluate(0.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N31, IsZeroAt1_0) {
  ASSERT_THAT(basis_function_.Evaluate(1.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N31, IsZeroAt1_5) {
  ASSERT_THAT(basis_function_.Evaluate(1.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N31, IsZeroAtMinus1_5) {
  ASSERT_THAT(basis_function_.Evaluate(-1.5), DoubleEq(0.0));
}

// Test basis function N_{0,2} from NURBS book example 2.1
class BasisFunctionEx21N02 : public Test {
 public:
  BasisFunctionEx21N02() :
      knot_vector_(vector<double>({0, 0, 0, 1, 1, 1})),
      basis_function_(knot_vector_, 2, 0) {}

 protected:
  KnotVector knot_vector_;
  BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N02, IsZeroAt0_0) {
  ASSERT_THAT(basis_function_.Evaluate(0.0), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx21N02, IsZeroAt0_5) {
  ASSERT_THAT(basis_function_.Evaluate(0.5), DoubleEq(0.25));
}

TEST_F(BasisFunctionEx21N02, IsZeroAt1_0) {
  ASSERT_THAT(basis_function_.Evaluate(1.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N02, IsZeroAt1_5) {
  ASSERT_THAT(basis_function_.Evaluate(1.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N02, IsZeroAtMinus1_5) {
  ASSERT_THAT(basis_function_.Evaluate(-1.5), DoubleEq(0.0));
}

// Test basis function N_{1,2} from NURBS book example 2.1
class BasisFunctionEx21N12 : public Test {
 public:
  BasisFunctionEx21N12() :
      knot_vector_(vector<double>({0, 0, 0, 1, 1, 1})),
      basis_function_(knot_vector_, 2, 1) {}

 protected:
  KnotVector knot_vector_;
  BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N12, IsZeroAt0_0) {
  ASSERT_THAT(basis_function_.Evaluate(0.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N12, IsZeroAt0_5) {
  ASSERT_THAT(basis_function_.Evaluate(0.5), DoubleEq(0.5));
}

TEST_F(BasisFunctionEx21N12, IsZeroAt1_0) {
  ASSERT_THAT(basis_function_.Evaluate(1.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N12, IsZeroAt1_5) {
  ASSERT_THAT(basis_function_.Evaluate(1.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N12, IsZeroAtMinus1_5) {
  ASSERT_THAT(basis_function_.Evaluate(-1.5), DoubleEq(0.0));
}

// Test basis function N_{2,2} from NURBS book example 2.1
class BasisFunctionEx21N22 : public Test {
 public:
  BasisFunctionEx21N22() :
      knot_vector_(vector<double>({0, 0, 0, 1, 1, 1})),
      basis_function_(knot_vector_, 2, 2) {}

 protected:
  KnotVector knot_vector_;
  BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N22, IsZeroAt0_0) {
  ASSERT_THAT(basis_function_.Evaluate(0.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N22, IsZeroAt0_5) {
  ASSERT_THAT(basis_function_.Evaluate(0.5), DoubleEq(0.25));
}

TEST_F(BasisFunctionEx21N22, IsOneAt1_0) {
  ASSERT_THAT(basis_function_.Evaluate(1.0), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx21N22, IsZeroAt1_5) {
  ASSERT_THAT(basis_function_.Evaluate(1.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N22, IsZeroAtMinus1_5) {
  ASSERT_THAT(basis_function_.Evaluate(-1.5), DoubleEq(0.0));
}

// Test basis function N_{0,1} from NURBS book example 2.2
class BasisFunctionEx22N01 : public Test {
 public:
  BasisFunctionEx22N01() : knot_vector_({0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}),
                           basis_function_(knot_vector_, 1, 0) {}

 protected:
  KnotVector knot_vector_;
  BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt0_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(0, 0), DoubleEq(basis_function_.Evaluate(0.0)));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt1_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(0, 1.5), DoubleEq(basis_function_.Evaluate(1.5)));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt2_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(0, 2.0), DoubleEq(basis_function_.Evaluate(2.0)));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt4_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(0, 4.0), DoubleEq(basis_function_.Evaluate(4.0)));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt5_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(0, 5.0), DoubleEq(basis_function_.Evaluate(5.0)));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt6_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(0, 6.0), DoubleEq(basis_function_.Evaluate(6.0)));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAtMinus0_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(0, -0.5),
              DoubleEq(basis_function_.Evaluate(-0.5)));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt0_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 0.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt1_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 1.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt2_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 2.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt4_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 4.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt5_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 5.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt6_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 6.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAtMinus0_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, -0.5), DoubleEq(0.0));
}

// Test basis function derivative N_{3,1} from NURBS book example 2.2
class BasisFunctionEx22N13 : public Test {
 public:
  BasisFunctionEx22N13() : knot_vector_({0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}),
                           basis_function_(knot_vector_, 1, 3) {}

 protected:
  KnotVector knot_vector_;
  BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqualZeroAt0_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 0.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqualZeroAt0_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 0.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqual1At1_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 1.0), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqual1At1_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 1.5), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqualMinus1At2_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 2.0), DoubleEq(-1.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqualMinus1At2_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 2.5), DoubleEq(-1.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqual0At3_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 3.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N13, FourthDerevitveIsEqual0At1_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(4, 1.5), DoubleEq(0.0));
}

// Test basis function derivative N_{6,1} from NURBS book example 2.2
class BasisFunctionEx22N61 : public Test {
 public:
  BasisFunctionEx22N61() : knot_vector_({0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}),
                           basis_function_(knot_vector_, 1, 6) {}

 protected:
  KnotVector knot_vector_;
  BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx22N61, FirstDerevitveIsEqualMinus1At4_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 4.0), DoubleEq(-1.0));
}

TEST_F(BasisFunctionEx22N61, SecondDerevitveIsEqual0At4_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(2, 4.5), DoubleEq(0.0));
}

// Test basis function derivative N_{7,2} from NURBS book example 2.2
class BasisFunctionEx22N72 : public Test {
 public:
  BasisFunctionEx22N72() : knot_vector_({0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5}),
                           basis_function_(knot_vector_, 2, 7) {}

 protected:
  KnotVector knot_vector_;
  BSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual0At4_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 4.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, SecondDerevitveIsEqual2At4_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(2, 4.0), DoubleEq(2.0));
}

TEST_F(BasisFunctionEx22N72, ThirdDerevitveIsEqual0At4_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(3, 4.0), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual1At4_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 4.5), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx22N72, SecondDerevitveIsEqual2At4_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(2, 4.5), DoubleEq(2.0));
}

TEST_F(BasisFunctionEx22N72, ThirdDerevitveIsEqual0At4_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(3, 4.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual2At5_0) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 5.0), DoubleEq(2.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual0At5_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 5.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual0At1_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, 1.5), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual0AtMinus5_5) {
  ASSERT_THAT(basis_function_.EvaluateDerivative(1, -5.5), DoubleEq(0.0));
}

