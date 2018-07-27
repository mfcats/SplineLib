/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "b_spline_basis_function.h"

#include "basis_function_factory.h"

baf::BSplineBasisFunction::BSplineBasisFunction(const KnotVector &knot_vector, Degree deg, uint64_t start_of_support)
    : BasisFunction(knot_vector, deg, KnotSpan{start_of_support}), start_knot_(knot_vector.GetKnot(start_of_support)), end_knot_(knot_vector.GetKnot(start_of_support+deg.get()+1)), left_denom_inv_(InverseWithPossiblyZeroDenominator(knot_vector.GetKnot(start_of_support+deg.get()).get() - start_knot_.get())), right_denom_inv_(InverseWithPossiblyZeroDenominator(end_knot_.get() - knot_vector.GetKnot(start_of_support+1).get())) {
  SetLowerDegreeBasisFunctions(knot_vector, start_of_support, deg);
}

double baf::BSplineBasisFunction::EvaluateOnSupport(ParamCoord param_coord) const {
  return ComputeLeftQuotient(param_coord) * left_lower_degree_->Evaluate(param_coord)
      + ComputeRightQuotient(param_coord) * right_lower_degree_->Evaluate(param_coord);
}

double baf::BSplineBasisFunction::EvaluateDerivativeOnSupport(ParamCoord param_coord, Derivative derivative) const {
  return GetDegree().get()
      * (ComputeLeftQuotientDenominatorInverse() * left_lower_degree_->EvaluateDerivative(param_coord, derivative - Derivative{1})
          - ComputeRightQuotientDenominatorInverse()
              * right_lower_degree_->EvaluateDerivative(param_coord, derivative - Derivative{1}));
}

void baf::BSplineBasisFunction::SetLowerDegreeBasisFunctions(const KnotVector &knot_vector,
                                                             uint64_t start_of_support,
                                                             Degree deg) {
  BasisFunctionFactory factory;

  left_lower_degree_.reset(factory.CreateDynamic(knot_vector, start_of_support, deg - Degree{1}));
  right_lower_degree_.reset(factory.CreateDynamic(knot_vector, start_of_support + 1, deg - Degree{1}));
}

double baf::BSplineBasisFunction::ComputeLeftQuotientDenominatorInverse() const {
  return left_denom_inv_;
}

double baf::BSplineBasisFunction::ComputeRightQuotientDenominatorInverse() const {
  return right_denom_inv_;
}

double baf::BSplineBasisFunction::ComputeLeftQuotient(ParamCoord param_coord) const {
  return (param_coord - start_knot_).get() * ComputeLeftQuotientDenominatorInverse();
}

double baf::BSplineBasisFunction::ComputeRightQuotient(ParamCoord param_coord) const {
  return (end_knot_ - param_coord).get() * ComputeRightQuotientDenominatorInverse();
}

double baf::BSplineBasisFunction::InverseWithPossiblyZeroDenominator(double denominator) {
  return std::abs(denominator) < util::NumericSettings<double>::kEpsilon() ? 0.0 : 1.0 / denominator;
}

