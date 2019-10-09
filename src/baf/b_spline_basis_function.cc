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
#include "src/util/numeric_settings.h"

namespace splinelib::src::baf {
BSplineBasisFunction::BSplineBasisFunction(const KnotVector &knot_vector,
                                                const Degree &degree,
                                                const KnotSpan &start_of_support)
    : BasisFunction(knot_vector, degree, start_of_support) {
  auto start_index = static_cast<size_t>(start_of_support.Get());
  auto degree_index = static_cast<size_t>(degree.Get());
  auto left_denom = (knot_vector.GetKnot(start_index + degree_index) - GetStartKnot()).Get();
  left_denom_inv_ = InverseWithPossiblyZeroDenominator(left_denom);
  auto right_denom = (GetEndKnot() - knot_vector.GetKnot(start_index + 1)).Get();
  right_denom_inv_ = InverseWithPossiblyZeroDenominator(right_denom);
  SetLowerDegreeBasisFunctions(knot_vector, degree, start_of_support);
}

double BSplineBasisFunction::EvaluateOnSupport(const ParametricCoordinate &param_coord) const {
  return ComputeLeftQuotient(param_coord) * left_lower_degree_->Evaluate(param_coord)
      + ComputeRightQuotient(param_coord) * right_lower_degree_->Evaluate(param_coord);
}

double BSplineBasisFunction::EvaluateDerivativeOnSupport(const ParametricCoordinate &param_coord,
                                                              const Derivative &derivative) const {
  return GetDegree().Get()
      * (left_denom_inv_ * left_lower_degree_->EvaluateDerivative(param_coord, derivative - Derivative{1})
          - right_denom_inv_ * right_lower_degree_->EvaluateDerivative(param_coord, derivative - Derivative{1}));
}

void BSplineBasisFunction::SetLowerDegreeBasisFunctions(const KnotVector &knot_vector,
                                                             const Degree &degree,
                                                             const KnotSpan &start_of_support) {
  left_lower_degree_.reset(BasisFunctionFactory::CreateDynamic(knot_vector, start_of_support, degree - Degree{1}));
  right_lower_degree_.reset(BasisFunctionFactory::CreateDynamic(knot_vector,
                                                                start_of_support + KnotSpan{1},
                                                                degree - Degree{1}));
}

double BSplineBasisFunction::ComputeLeftQuotient(const ParametricCoordinate &param_coord) const {
  return (param_coord - GetStartKnot()).Get() * left_denom_inv_;
}

double BSplineBasisFunction::ComputeRightQuotient(const ParametricCoordinate &param_coord) const {
  return (GetEndKnot() - param_coord).Get() * right_denom_inv_;
}

double BSplineBasisFunction::InverseWithPossiblyZeroDenominator(double denominator) const {
  return std::abs(denominator) < util::numeric_settings::GetEpsilon<double>() ? 0.0 : 1.0 / denominator;
}
}  // namespace splinelib::src::baf
