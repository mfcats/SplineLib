/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_BAF_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_H_
#define SRC_BAF_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_H_

#include "src/baf/b_spline_basis_function.h"
#include "src/baf/knot_vector.h"
#include "src/util/named_type.h"

// A ZeroDegreeBSplineBasisFunction N_{i,0} is a step function. It equals one in the half-open interval [u_i, u_{i+p+1})
// which is stored in the member variables start_knot_ and end_knot_ of the base class BSplineBasisFunction resp. in the
// closed interval [u_i, u_{i+p+1}] when end_knot_ u_{i+1} equals the last knot u_m of the knot vector and then
// end_knot_is_last_knot_ is set to true.
// Everywhere else the function equals zero.
// The derivative is defined to be zero on the whole knot vector interval [u_0, u_m].
// Example (basis function N_{2,0} for simplest knot vector of degree 2, see NURBS book example 2.1):
//   ZeroDegreeBSplineBasisFunction basis_function(KnotVector{0.0, 0.0, 0.0, 1.0, 1.0, 1.0}, KnotSpan{2});
//   Degree p = basis_function.GetDegree();  // Returns 0.
//   ParametricCoordinate start_knot = basis_function.GetStartKnot();  // Returns 0.0;
//   ParametricCoordinate end_knot = basis_function.GetEndKnot();  // Returns 1.0;
//   double evaluate = basis_function.Evaluate(ParametricCoordinate{0.5});  // Returns 1.0;
//   double evaluate_derivative = basis_function.EvaluateDerivative(ParametricCoordinate{0.5});  // Returns 0.0;
//   // The basis functions with knot span 0, 1 or 3 are always evaluated to zero as they are defined on an interval of
//   // length zero. For N_{3,0} end_knot_is_last_knot_ would equal true.
namespace splinelib::src::baf {
class ZeroDegreeBSplineBasisFunction : public BSplineBasisFunction {
 public:
  ZeroDegreeBSplineBasisFunction() = default;
  ZeroDegreeBSplineBasisFunction(KnotVector const &knot_vector, KnotSpan const &start_of_support);
  ZeroDegreeBSplineBasisFunction(ZeroDegreeBSplineBasisFunction const &other) = delete;
  ZeroDegreeBSplineBasisFunction(ZeroDegreeBSplineBasisFunction &&other) = delete;
  ZeroDegreeBSplineBasisFunction & operator=(ZeroDegreeBSplineBasisFunction const &rhs) = delete;
  ZeroDegreeBSplineBasisFunction & operator=(ZeroDegreeBSplineBasisFunction &&rhs) = delete;
  ~ZeroDegreeBSplineBasisFunction() final = default;

 protected:
  double EvaluateOnSupport(ParametricCoordinate const &/*parametric_coordinate*/) const final;
  double EvaluateDerivativeOnSupport(ParametricCoordinate const &/*parametric_coordinate*/,
                                     Derivative const &/*derivative*/) const final;
};

#include "src/baf/zero_degree_b_spline_basis_function.inc"
}  // namespace splinelib::src::baf

#endif  // SRC_BAF_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_H_
