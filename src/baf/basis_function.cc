/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "basis_function.h"

double baf::BasisFunction::Evaluate(const ParamCoord &paramCoord) const {
  return IsCoordinateInSupport(paramCoord) ? this->EvaluateOnSupport(paramCoord) : 0.0;
}

double baf::BasisFunction::EvaluateDerivative(const ParamCoord &param_coord, const Derivative &derivative) const {
  return derivative.get() == 0 ? Evaluate(param_coord) :
         (IsCoordinateInSupport(param_coord) ? this->EvaluateDerivativeOnSupport(param_coord, derivative) : 0.0);
}

baf::BasisFunction::BasisFunction(const KnotVector &knot_vector, const Degree &degree, const KnotSpan &start_of_support)
    : degree_(degree) {
  auto start_index = static_cast<size_t>(start_of_support.get());
  auto degree_index = static_cast<size_t>(degree.get());
  start_knot_ = knot_vector.GetKnot(start_index);
  end_knot_ = knot_vector.GetKnot(start_index + degree_index + 1);
  end_knot_is_last_knot_ = knot_vector.IsLastKnot(end_knot_);
}

Degree baf::BasisFunction::GetDegree() const {
  return degree_;
}

ParamCoord baf::BasisFunction::GetStartKnot() const {
  return start_knot_;
}

ParamCoord baf::BasisFunction::GetEndKnot() const {
  return end_knot_;
}

bool baf::BasisFunction::IsCoordinateInSupport(const ParamCoord &param_coord) const {
  return (start_knot_ <= param_coord && param_coord < end_knot_)
      || (end_knot_is_last_knot_ && param_coord == end_knot_);
}
