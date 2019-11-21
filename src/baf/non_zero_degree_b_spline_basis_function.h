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

namespace splinelib::src::baf {
class NonZeroDegreeBSplineBasisFunction : public BSplineBasisFunction {
 public:
  NonZeroDegreeBSplineBasisFunction(const KnotVector &knot_vector, const Degree &degree,
                                    const KnotSpan &start_of_support);

 protected:
  double EvaluateOnSupport(const ParametricCoordinate &param_coord) const override;
  double EvaluateDerivativeOnSupport(const ParametricCoordinate &param_coord,
                                     const Derivative &derivative) const override;

 private:
  void SetLowerDegreeBasisFunctions(const KnotVector &knot_vector,
                                    const Degree &degree,
                                    const KnotSpan &start_of_support);

  double ComputeLeftQuotient(const ParametricCoordinate &param_coord) const;
  double ComputeRightQuotient(const ParametricCoordinate &param_coord) const;

  double InverseWithPossiblyZeroDenominator(double denominator) const;

  double left_denom_inv_;
  double right_denom_inv_;
  std::unique_ptr<BSplineBasisFunction> left_lower_degree_;
  std::unique_ptr<BSplineBasisFunction> right_lower_degree_;
};
}  // namespace splinelib::src::baf

#endif  // SRC_BAF_NON_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_H_
