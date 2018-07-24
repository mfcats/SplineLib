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

#include "numeric_settings.h"

double baf::BasisFunction::Evaluate(ParamCoord paramCoord) const {
  return IsCoordinateInSupport(paramCoord) ? this->EvaluateOnSupport(paramCoord) : 0.0;
}

double baf::BasisFunction::EvaluateDerivative(ParamCoord param_coord, int derivative) const {
  return derivative == 0 ? Evaluate(param_coord) :
         IsCoordinateInSupport(param_coord) ? this->EvaluateDerivativeOnSupport(param_coord, derivative) : 0.0;
}

baf::BasisFunction::BasisFunction(KnotVector knot_vector, int degree, uint64_t start)
    : knotVector_(std::move(knot_vector)), degree_(degree), start_of_support_(start) {}

ParamCoord baf::BasisFunction::GetKnot(uint64_t knot_position) const {
  return knotVector_.GetKnot(knot_position);
}

uint64_t baf::BasisFunction::GetStartOfSupport() const {
  return start_of_support_;
}

int baf::BasisFunction::GetDegree() const {
  return degree_;
}

bool baf::BasisFunction::IsCoordinateInSupport(ParamCoord param_coord) const {
  return knotVector_.IsInKnotVectorRange(param_coord) && IsCoordinateInSupportSpan(param_coord);
}

bool baf::BasisFunction::IsCoordinateInSupportSpan(ParamCoord param_coord) const {
  auto span = knotVector_.GetKnotSpan(param_coord);
  return !(span < start_of_support_ || span >= start_of_support_ + degree_ + 1);
}
