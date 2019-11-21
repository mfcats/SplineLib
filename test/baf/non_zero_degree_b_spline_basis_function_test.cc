/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/baf/non_zero_degree_b_spline_basis_function.h"
#include "src/util/stl_container_access.h"

using testing::DoubleEq;
using testing::Test;

using namespace splinelib::src;

class MockKnotVector000111 : public baf::KnotVector {
 public:
  ParametricCoordinate operator[](int knot_num) const override {
    return GetValue(knots_, knot_num);
  }

  bool IsLastKnot(const ParametricCoordinate &param_coord) const override {
    return param_coord == ParametricCoordinate{1};
  }

 private:
  const std::vector<ParametricCoordinate> knots_ = {ParametricCoordinate{0}, ParametricCoordinate{0},
                                                    ParametricCoordinate{0}, ParametricCoordinate{1},
                                                    ParametricCoordinate{1}, ParametricCoordinate{1}};
};

// Test basis function N_{0,1} from NURBS book example 2.1
class BasisFunctionEx21N01 : public Test {
 public:
  BasisFunctionEx21N01() : basis_function_(knot_vector_000111, Degree{1}, KnotSpan{0}) {}

 protected:
  const MockKnotVector000111 knot_vector_000111;
  baf::NonZeroDegreeBSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N01, IsZeroAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N01, IsZeroAt0_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N01, IsZeroAt1_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N01, IsZeroAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N01, IsZeroAtMinus1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{1,1} from NURBS book example 2.1
class BasisFunctionEx21N11 : public Test {
 public:
  BasisFunctionEx21N11() : basis_function_(knot_vector_000111, Degree{1}, KnotSpan{1}) {}

 protected:
  const MockKnotVector000111 knot_vector_000111;
  baf::NonZeroDegreeBSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N11, IsOneAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.0}), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx21N11, Is0_5At0_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.5}), DoubleEq(0.5));
}

TEST_F(BasisFunctionEx21N11, IsZeroAt1_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N11, IsZeroAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N11, IsZeroAtMinus1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{2,1} from NURBS book example 2.1
class BasisFunctionEx21N21 : public Test {
 public:
  BasisFunctionEx21N21() : basis_function_(knot_vector_000111, Degree{1}, KnotSpan{2}) {}

 protected:
  const MockKnotVector000111 knot_vector_000111;
  baf::NonZeroDegreeBSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N21, IsZeroAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N21, Is0_5At0_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.5}), DoubleEq(0.5));
}

TEST_F(BasisFunctionEx21N21, IsOneAt1_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.0}), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx21N21, IsZeroAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N21, IsZeroAtMinus1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{3,1} from NURBS book example 2.1
class BasisFunctionEx21N31 : public Test {
 public:
  BasisFunctionEx21N31() : basis_function_(knot_vector_000111, Degree{1}, KnotSpan{3}) {}

 protected:
  const MockKnotVector000111 knot_vector_000111;
  baf::NonZeroDegreeBSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N31, IsZeroAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N31, IsZeroAt0_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N31, IsZeroAt1_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N31, IsZeroAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N31, IsZeroAtMinus1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{0,2} from NURBS book example 2.1
class BasisFunctionEx21N02 : public Test {
 public:
  BasisFunctionEx21N02() : basis_function_(knot_vector_000111, Degree{2}, KnotSpan{0}) {}

 protected:
  const MockKnotVector000111 knot_vector_000111;
  baf::NonZeroDegreeBSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N02, IsZeroAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.0}), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx21N02, Is0_25At0_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.5}), DoubleEq(0.25));
}

TEST_F(BasisFunctionEx21N02, IsZeroAt1_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N02, IsZeroAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N02, IsZeroAtMinus1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{1,2} from NURBS book example 2.1
class BasisFunctionEx21N12 : public Test {
 public:
  BasisFunctionEx21N12() : basis_function_(knot_vector_000111, Degree{2}, KnotSpan{1}) {}

 protected:
  const MockKnotVector000111 knot_vector_000111;
  baf::NonZeroDegreeBSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N12, IsZeroAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N12, Is0_5At0_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.5}), DoubleEq(0.5));
}

TEST_F(BasisFunctionEx21N12, IsZeroAt1_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N12, IsZeroAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N12, IsZeroAtMinus1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{2,2} from NURBS book example 2.1
class BasisFunctionEx21N22 : public Test {
 public:
  BasisFunctionEx21N22() : basis_function_(knot_vector_000111, Degree{2}, KnotSpan{2}) {}

 protected:
  const MockKnotVector000111 knot_vector_000111;
  baf::NonZeroDegreeBSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx21N22, IsZeroAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.0}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N22, Is0_25At0_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{0.5}), DoubleEq(0.25));
}

TEST_F(BasisFunctionEx21N22, IsOneAt1_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.0}), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx21N22, IsZeroAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{1.5}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx21N22, IsZeroAtMinus1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParametricCoordinate{-1.5}), DoubleEq(0.0));
}

class MockKnotVector00012344555 : public baf::KnotVector {
 public:
  ParametricCoordinate operator[](int knot_num) const override {
    return GetValue(knots_, knot_num);
  }

  bool IsLastKnot(const ParametricCoordinate &param_coord) const override {
    return param_coord == ParametricCoordinate{5};
  }

 private:
  const std::vector<ParametricCoordinate> knots_ =
      {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1},
       ParametricCoordinate{2}, ParametricCoordinate{3}, ParametricCoordinate{4}, ParametricCoordinate{4},
       ParametricCoordinate{5}, ParametricCoordinate{5}, ParametricCoordinate{5}};
};

// Test basis function N_{0,1} from NURBS book example 2.2
class BasisFunctionEx22N01 : public Test {
 public:
  BasisFunctionEx22N01() : basis_function_(knot_vector_00012344555, Degree{1}, KnotSpan{0}) {}

 protected:
  const MockKnotVector00012344555 knot_vector_00012344555;
  baf::NonZeroDegreeBSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{0}, Derivative{0}),
              DoubleEq(basis_function_.Evaluate(ParametricCoordinate{0.0})));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{1.5}, Derivative{0}),
              DoubleEq(basis_function_.Evaluate(ParametricCoordinate{1.5})));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt2_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{2.0}, Derivative{0}),
              DoubleEq(basis_function_.Evaluate(ParametricCoordinate{2.0})));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt4_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{4.0}, Derivative{0}),
              DoubleEq(basis_function_.Evaluate(ParametricCoordinate{4.0})));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt5_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{5.0}, Derivative{0}),
              DoubleEq(basis_function_.Evaluate(ParametricCoordinate{5.0})));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAt6_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{6.0}, Derivative{0}),
              DoubleEq(basis_function_.Evaluate(ParametricCoordinate{6.0})));
}

TEST_F(BasisFunctionEx22N01, ZerothDerevitveIsEqualValueAtMinus0_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{-0.5}, Derivative{0}),
              DoubleEq(basis_function_.Evaluate(ParametricCoordinate{-0.5})));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{0.0}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{1.5}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt2_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{2.0}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt4_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{4.0}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt5_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{5.0}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAt6_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{6.0}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N01, FirstDerevitveIsEqualZeroAtMinus0_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{-0.5}, Derivative{1}), DoubleEq(0.0));
}

// Test basis function derivative N_{3,1} from NURBS book example 2.2
class BasisFunctionEx22N13 : public Test {
 public:
  BasisFunctionEx22N13() : basis_function_(knot_vector_00012344555, Degree{1}, KnotSpan{3}) {}

 protected:
  const MockKnotVector00012344555 knot_vector_00012344555;
  baf::NonZeroDegreeBSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqualZeroAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{0.0}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqualZeroAt0_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{0.5}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqual1At1_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{1.0}, Derivative{1}), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqual1At1_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{1.5}, Derivative{1}), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqualMinus1At2_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{2.0}, Derivative{1}), DoubleEq(-1.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqualMinus1At2_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{2.5}, Derivative{1}), DoubleEq(-1.0));
}

TEST_F(BasisFunctionEx22N13, FirstDerevitveIsEqual0At3_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{3.0}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N13, FourthDerevitveIsEqual0At1_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{1.5}, Derivative{4}), DoubleEq(0.0));
}

// Test basis function derivative N_{6,1} from NURBS book example 2.2
class BasisFunctionEx22N61 : public Test {
 public:
  BasisFunctionEx22N61() : basis_function_(knot_vector_00012344555, Degree{1}, KnotSpan{6}) {}

 protected:
  const MockKnotVector00012344555 knot_vector_00012344555;
  baf::NonZeroDegreeBSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx22N61, FirstDerevitveIsEqualMinus1At4_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{4.0}, Derivative{1}), DoubleEq(-1.0));
}

TEST_F(BasisFunctionEx22N61, SecondDerevitveIsEqual0At4_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{4.5}, Derivative{2}), DoubleEq(0.0));
}

// Test basis function derivative N_{7,2} from NURBS book example 2.2
class BasisFunctionEx22N72 : public Test {
 public:
  BasisFunctionEx22N72() : basis_function_(knot_vector_00012344555, Degree{2}, KnotSpan{7}) {}

 protected:
  const MockKnotVector00012344555 knot_vector_00012344555;
  baf::NonZeroDegreeBSplineBasisFunction basis_function_;
};

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual0At4_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{4.0}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, SecondDerevitveIsEqual2At4_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{4.0}, Derivative{2}), DoubleEq(2.0));
}

TEST_F(BasisFunctionEx22N72, ThirdDerevitveIsEqual0At4_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{4.0}, Derivative{3}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual1At4_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{4.5}, Derivative{1}), DoubleEq(1.0));
}

TEST_F(BasisFunctionEx22N72, SecondDerevitveIsEqual2At4_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{4.5}, Derivative{2}), DoubleEq(2.0));
}

TEST_F(BasisFunctionEx22N72, ThirdDerevitveIsEqual0At4_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{4.5}, Derivative{3}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual2At5_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{5.0}, Derivative{1}), DoubleEq(2.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual0At5_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{5.5}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual0At1_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{1.5}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(BasisFunctionEx22N72, FirstDerevitveIsEqual0AtMinus5_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParametricCoordinate{-5.5}, Derivative{1}), DoubleEq(0.0));
}
