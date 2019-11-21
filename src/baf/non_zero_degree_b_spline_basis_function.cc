/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "src/baf/non_zero_degree_b_spline_basis_function.h"

#include "src/util/numeric_settings.h"

namespace splinelib::src::baf {
NonZeroDegreeBSplineBasisFunction::NonZeroDegreeBSplineBasisFunction(KnotVector const &knot_vector,
                                                                     KnotSpan const &start_of_support,
                                                                     Degree const &degree)
    : BSplineBasisFunction(knot_vector, start_of_support, degree) {
  auto const start_index = start_of_support.Get();
  auto const degree_index = degree.Get();
  auto const left_denominator = (knot_vector[start_index + degree_index] - GetStartKnot()).Get();
  left_denominator_inverse_ = InverseWithPossiblyZeroDenominator(left_denominator);
  auto const right_denominator = (BSplineBasisFunction::GetEndKnot() - knot_vector[start_index + 1]).Get();
  right_denominator_inverse_ = InverseWithPossiblyZeroDenominator(right_denominator);
  SetLowerDegreeBasisFunctions(knot_vector, start_of_support, degree);
}

double NonZeroDegreeBSplineBasisFunction::EvaluateOnSupport(ParametricCoordinate const &parametric_coordinate) const {
  return (ComputeLeftQuotient(parametric_coordinate) *
          left_lower_degree_basis_function_->Evaluate(parametric_coordinate) +
          ComputeRightQuotient(parametric_coordinate) *
          right_lower_degree_basis_function_->Evaluate(parametric_coordinate));
}

double NonZeroDegreeBSplineBasisFunction::EvaluateDerivativeOnSupport(ParametricCoordinate const &parametric_coordinate,
                                                                      Derivative const &derivative) const {
  Derivative const derivative_reduced_by_one = derivative - Derivative{1};
  return (BSplineBasisFunction::GetDegree().Get() *
          (left_denominator_inverse_ *
           left_lower_degree_basis_function_->EvaluateDerivative(parametric_coordinate, derivative_reduced_by_one) -
           right_denominator_inverse_ *
           right_lower_degree_basis_function_->EvaluateDerivative(parametric_coordinate, derivative_reduced_by_one)));
}

void NonZeroDegreeBSplineBasisFunction::SetLowerDegreeBasisFunctions(KnotVector const &knot_vector,
                                                                     KnotSpan const &start_of_support,
                                                                     Degree const &degree) {
  Degree const degree_reduced_by_one = degree - Degree{1};
  left_lower_degree_basis_function_.reset(
      BSplineBasisFunction::CreateDynamic(knot_vector, start_of_support, degree_reduced_by_one));
  right_lower_degree_basis_function_.reset(
      BSplineBasisFunction::CreateDynamic(knot_vector, start_of_support + KnotSpan{1}, degree_reduced_by_one));
}

double NonZeroDegreeBSplineBasisFunction::InverseWithPossiblyZeroDenominator(double denominator) const {
  if (std::abs(denominator) < util::numeric_settings::GetEpsilon<double>()) return 0.0;
  return (1.0 / denominator);
}
}  // namespace splinelib::src::baf
