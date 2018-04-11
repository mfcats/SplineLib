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

BSplineBasisFunction::BSplineBasisFunction(KnotVector knot_vector,
                                           Degree deg,
                                           uint64_t start_of_support) :
    BasisFunction(knot_vector, deg, start_of_support) {
  SetLowerDegreeBsisFunctions(knot_vector, start_of_support, deg);
}

double BSplineBasisFunction::EvaluateOnSupport(double param_coord) const {
  return ComputeLeftQuotient(param_coord) * left_lower_degree_->Evaluate(param_coord) +
      ComputeRightQuotient(param_coord) * right_lower_degree_->Evaluate(param_coord);
}

double BSplineBasisFunction::EvaluateDerivativeOnSupport(Derivative derivative,
                                                         double param_coord) const {
  return GetDegree() * (ComputeLeftQuotientDenominatorInverse()
      * left_lower_degree_->EvaluateDerivative(derivative - 1, param_coord) -
      ComputeRightQuotientDenominatorInverse()
          * right_lower_degree_->EvaluateDerivative(derivative - 1, param_coord));
}

void BSplineBasisFunction::SetLowerDegreeBsisFunctions(KnotVector knot_vector,
                                                       uint64_t start_of_support,
                                                       Degree deg) {
  BasisFunctionFactory basis_function_factory;

  left_lower_degree_.reset(basis_function_factory.CreateDynamic(knot_vector,
                                                                start_of_support,
                                                                deg - 1));
  right_lower_degree_.reset(basis_function_factory.CreateDynamic(knot_vector,
                                                                 start_of_support + 1,
                                                                 deg - 1));
}

double BSplineBasisFunction::ComputeLeftQuotientDenominatorInverse() const {
  auto left_quotient_denominator =
      GetKnot(GetStartOfSupport() + GetDegree()) - GetKnot(GetStartOfSupport());
  return std::fabs(left_quotient_denominator) < NumericSettings<double>::kEpsilon() ? 0.0 : 1.0
      / left_quotient_denominator;
}

double BSplineBasisFunction::ComputeRightQuotientDenominatorInverse() const {
  auto right_quotient_denominator =
      GetKnot(GetStartOfSupport() + GetDegree() + 1) - GetKnot(GetStartOfSupport() + 1);
  return std::fabs(right_quotient_denominator) < NumericSettings<double>::kEpsilon() ? 0.0 : 1.0
      / right_quotient_denominator;
}

double BSplineBasisFunction::ComputeLeftQuotient(double param_coord) const {
  return (param_coord - GetKnot(GetStartOfSupport())) * ComputeLeftQuotientDenominatorInverse();
}
double BSplineBasisFunction::ComputeRightQuotient(double param_coord) const {
  return (GetKnot(GetStartOfSupport() + GetDegree() + 1) - param_coord)
      * ComputeRightQuotientDenominatorInverse();
}

