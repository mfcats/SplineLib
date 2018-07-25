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

baf::BSplineBasisFunction::BSplineBasisFunction(const KnotVector &knot_vector, int deg, uint64_t start_of_support)
    : BasisFunction(knot_vector, deg, start_of_support) {
  SetLowerDegreeBasisFunctions(knot_vector, start_of_support, deg);
}

double baf::BSplineBasisFunction::EvaluateOnSupport(ParamCoord param_coord) const {
  return ComputeLeftQuotient(param_coord) * left_lower_degree_->Evaluate(param_coord, knotVector_)
      + ComputeRightQuotient(param_coord) * right_lower_degree_->Evaluate(param_coord, knotVector_);
}

double baf::BSplineBasisFunction::EvaluateDerivativeOnSupport(ParamCoord param_coord, int derivative) const {
  return GetDegree().get()
      * (ComputeLeftQuotientDenominatorInverse() * left_lower_degree_->EvaluateDerivative(param_coord, derivative - 1)
          - ComputeRightQuotientDenominatorInverse()
              * right_lower_degree_->EvaluateDerivative(param_coord, derivative - 1));
}

void baf::BSplineBasisFunction::SetLowerDegreeBasisFunctions(const KnotVector &knot_vector,
                                                             uint64_t start_of_support,
                                                             int deg) {
  BasisFunctionFactory factory;

  left_lower_degree_.reset(factory.CreateDynamic(knot_vector, start_of_support, deg - 1));
  right_lower_degree_.reset(factory.CreateDynamic(knot_vector, start_of_support + 1, deg - 1));
}

double baf::BSplineBasisFunction::ComputeLeftQuotientDenominatorInverse() const {
  auto left_denom = GetKnot(GetStartOfSupport().get() + GetDegree().get()) - GetKnot(GetStartOfSupport().get());
  return InverseWithPossiblyZeroDenominator(left_denom.get());
}

double baf::BSplineBasisFunction::ComputeRightQuotientDenominatorInverse() const {
  auto right_denom = GetKnot(GetStartOfSupport().get() + GetDegree().get() + 1) - GetKnot(GetStartOfSupport().get() + 1);
  return InverseWithPossiblyZeroDenominator(right_denom.get());
}

double baf::BSplineBasisFunction::ComputeLeftQuotient(ParamCoord param_coord) const {
  return (param_coord - GetKnot(GetStartOfSupport().get())).get() * ComputeLeftQuotientDenominatorInverse();
}

double baf::BSplineBasisFunction::ComputeRightQuotient(ParamCoord param_coord) const {
  return (GetKnot(GetStartOfSupport().get() + GetDegree().get() + 1) - param_coord).get()
      * ComputeRightQuotientDenominatorInverse();
}

double baf::BSplineBasisFunction::InverseWithPossiblyZeroDenominator(double denominator) const {
  return std::fabs(denominator) < util::NumericSettings<double>::kEpsilon() ? 0.0 : 1.0 / denominator;
}

