/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_BAF_B_SPLINE_BASIS_FUNCTION_H_
#define SRC_BAF_B_SPLINE_BASIS_FUNCTION_H_

#include <cmath>
#include <limits>
#include <memory>

#include "basis_function.h"
#include "numeric_settings.h"

namespace baf {
class BSplineBasisFunction : public BasisFunction {
 public:
  BSplineBasisFunction(const KnotVector &knot_vector, Degree deg, uint64_t start_of_support);

 protected:
  double EvaluateOnSupport(ParamCoord param_coord) const override;

  double EvaluateDerivativeOnSupport(ParamCoord param_coord, Derivative derivative) const override;

  std::unique_ptr<BasisFunction> left_lower_degree_;
  std::unique_ptr<BasisFunction> right_lower_degree_;

 private:
  void SetLowerDegreeBasisFunctions(const KnotVector &knot_vector, uint64_t start_of_support, Degree deg);

  double ComputeLeftQuotientDenominatorInverse() const;

  double ComputeRightQuotientDenominatorInverse() const;

  double ComputeLeftQuotient(ParamCoord param_coord) const;

  double ComputeRightQuotient(ParamCoord param_coord) const;

  static double InverseWithPossiblyZeroDenominator(double denominator);

  ParamCoord start_knot_;
  ParamCoord end_knot_;
  double left_denom_inv_;
  double right_denom_inv_;

};
}  // namespace baf

#endif  // SRC_BAF_B_SPLINE_BASIS_FUNCTION_H_
