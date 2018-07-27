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

#include "zero_degree_b_spline_basis_function.h"

using std::vector;

using testing::DoubleEq;
using testing::Test;

// Test basis function N_{0,0} from NURBS book example 2.1
class ZeroDegreeBasisFunctionEx21N00 : public Test {
 public:
  ZeroDegreeBasisFunctionEx21N00() : knot_vector_(vector<ParamCoord>({ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                                      ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
                                     basis_function_(knot_vector_, KnotSpan{0}) {}

 protected:
  baf::KnotVector knot_vector_;
  baf::ZeroDegBSplBasFnc basis_function_;
};

TEST_F(ZeroDegreeBasisFunctionEx21N00, IsZeroAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.0}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx21N00, IsZeroAt0_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.5}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx21N00, IsZeroAt1_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.0}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx21N00, IsZeroAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.5}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx21N00, IsZeroAMinust1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{2,0} from NURBS book example 2.1
class ZeroDegreeBasisFunctionEx21N20 : public Test {
 public:
  ZeroDegreeBasisFunctionEx21N20() : knot_vector_(vector<ParamCoord>({ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                                      ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
                                     basis_function_(knot_vector_, KnotSpan{2}) {}

 protected:
  baf::KnotVector knot_vector_;
  baf::ZeroDegBSplBasFnc basis_function_;
};

TEST_F(ZeroDegreeBasisFunctionEx21N20, IsOneAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.0}), DoubleEq(1.0));
}

TEST_F(ZeroDegreeBasisFunctionEx21N20, IsOneAt0_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.5}), DoubleEq(1.0));
}

TEST_F(ZeroDegreeBasisFunctionEx21N20, IsOneAt1_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.0}), DoubleEq(1.0));
}

TEST_F(ZeroDegreeBasisFunctionEx21N20, IsZeroAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.5}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx21N20, IsZeroAMinust1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{4,0} from NURBS book example 2.1
class ZeroDegreeBasisFunctionEx21N40 : public Test {
 public:
  ZeroDegreeBasisFunctionEx21N40() : knot_vector_(vector<ParamCoord>({ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                                      ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
                                     basis_function_(knot_vector_, KnotSpan{0}) {}

 protected:
  baf::KnotVector knot_vector_;
  baf::ZeroDegBSplBasFnc basis_function_;
};

TEST_F(ZeroDegreeBasisFunctionEx21N40, IsZeroAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.0}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx21N40, IsZeroAt0_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{0.5}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx21N40, IsZeroAt1_0) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.0}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx21N40, IsZeroAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{1.5}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx21N40, IsZeroAMinust1_5) { // NOLINT
  ASSERT_THAT(basis_function_.Evaluate(ParamCoord{-1.5}), DoubleEq(0.0));
}

// Test basis function N_{0,0} from NURBS book example 2.2
class ZeroDegreeBasisFunctionEx22N00 : public Test {
 public:
  ZeroDegreeBasisFunctionEx22N00() : knot_vector_(vector<ParamCoord>({ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                                      ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
                                     basis_function_(knot_vector_, KnotSpan{0}) {}

 protected:
  baf::KnotVector knot_vector_;
  baf::ZeroDegBSplBasFnc basis_function_;
};

TEST_F(ZeroDegreeBasisFunctionEx22N00, ZerothDerevitveIsEqualValueAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{0}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{0.0})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N00, ZerothDerevitveIsEqualValueAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{1}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{1.5})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N00, ZerothDerevitveIsEqualValueAt2_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{2}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{2.0})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N00, ZerothDerevitveIsEqualValueAt4_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{4}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{4.0})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N00, ZerothDerevitveIsEqualValueAt5_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{5}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{5.0})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N00, ZerothDerevitveIsEqualValueAt6_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{6}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{6.0})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N00, ZerothDerevitveIsEqualValueAtMinus0_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{0}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{-0.5})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N00, FirstDerevitveIsEqualZeroAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{0}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx22N00, FirstDerevitveIsEqualZeroAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx22N00, FirstDerevitveIsEqualZeroAt2_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{2}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx22N00, FirstDerevitveIsEqualZeroAt4_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{4}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx22N00, FirstDerevitveIsEqualZeroAt5_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{5}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx22N00, FirstDerevitveIsEqualZeroAt6_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{6}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx22N00, FirstDerevitveIsEqualZeroAtMinus0_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{0}), DoubleEq(0.0));
}

// Test basis function N_{4,0} from NURBS book example 2.2
class ZeroDegreeBasisFunctionEx22N40 : public Test {
 public:
  ZeroDegreeBasisFunctionEx22N40() : knot_vector_(vector<ParamCoord>({ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                                      ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
                                     basis_function_(knot_vector_, KnotSpan{0}) {}

 protected:
  baf::KnotVector knot_vector_;
  baf::ZeroDegBSplBasFnc basis_function_;
};

TEST_F(ZeroDegreeBasisFunctionEx22N40, ZerothDerevitveIsEqualValueAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{0}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{0.0})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N40, ZerothDerevitveIsEqualValueAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{1}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{1.5})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N40, ZerothDerevitveIsEqualValueAt2_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{2}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{2.0})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N40, ZerothDerevitveIsEqualValueAt4_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{4}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{4.0})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N40, ZerothDerevitveIsEqualValueAt5_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{5}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{5.0})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N40, ZerothDerevitveIsEqualValueAt6_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{6}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{6.0})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N40, ZerothDerevitveIsEqualValueAtMinus0_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{0}, Derivative{0}),
              DoubleEq(basis_function_.Evaluate(ParamCoord{-0.5})));
}

TEST_F(ZeroDegreeBasisFunctionEx22N40, FirstDerevitveIsEqualZeroAt0_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{0}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx22N40, FirstDerevitveIsEqualZeroAt1_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{1}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx22N40, FirstDerevitveIsEqualZeroAt2_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{2}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx22N40, FirstDerevitveIsEqualZeroAt4_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{4}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx22N40, FirstDerevitveIsEqualZeroAt5_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{5}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx22N40, FirstDerevitveIsEqualZeroAt6_0) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{6}), DoubleEq(0.0));
}

TEST_F(ZeroDegreeBasisFunctionEx22N40, FirstDerevitveIsEqualZeroAtMinus0_5) { // NOLINT
  ASSERT_THAT(basis_function_.EvaluateDerivative(ParamCoord{1}, Derivative{0}), DoubleEq(0.0));
}
