/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_BAF_NON_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_H_
#define SRC_BAF_NON_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_H_

#include <memory>

#include "src/baf/b_spline_basis_function.h"

#include "src/baf/knot_vector.h"
#include "src/util/named_type.h"

namespace splinelib::src::baf {
class NonZeroDegreeBSplineBasisFunction : public BSplineBasisFunction {
 public:
  NonZeroDegreeBSplineBasisFunction(KnotVector const &knot_vector, KnotSpan const &start_of_support,
                                    Degree const &degree);
  NonZeroDegreeBSplineBasisFunction(NonZeroDegreeBSplineBasisFunction const &other) = delete;
  NonZeroDegreeBSplineBasisFunction(NonZeroDegreeBSplineBasisFunction &&other) = delete;
  NonZeroDegreeBSplineBasisFunction & operator=(NonZeroDegreeBSplineBasisFunction const &other) = delete;
  NonZeroDegreeBSplineBasisFunction & operator=(NonZeroDegreeBSplineBasisFunction &&other) = delete;
  ~NonZeroDegreeBSplineBasisFunction() final = default;

 protected:
  double EvaluateOnSupport(ParametricCoordinate const &parametric_coordinate) const final;
  double EvaluateDerivativeOnSupport(ParametricCoordinate const &parametric_coordinate,
                                     Derivative const &derivative) const final;

 private:
  void SetLowerDegreeBasisFunctions(KnotVector const &knot_vector, KnotSpan const &start_of_support,
                                    Degree const &degree);
  double ComputeLeftQuotient(ParametricCoordinate const &parametric_coordinate) const;
  double ComputeRightQuotient(ParametricCoordinate const &parametric_coordinate) const;

  // Within the context of the Cox-de Boor recursion formula, the inverse of 0 is defined to be 0 (see NURBS book
  // equation 2.5).
  double InverseWithPossiblyZeroDenominator(double denominator) const;

  double left_denominator_inverse_;
  double right_denominator_inverse_;
  std::unique_ptr<BSplineBasisFunction> left_lower_degree_basis_function_;
  std::unique_ptr<BSplineBasisFunction> right_lower_degree_basis_function_;
};

#include "src/baf/non_zero_degree_b_spline_basis_function.inc"
}  // namespace splinelib::src::baf

#endif  // SRC_BAF_NON_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_H_
