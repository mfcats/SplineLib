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

#include <cmath>

#include "numeric_settings.h"

double BasisFunction::Evaluate(double paramCoord) const {
  return IsCoordinateInSupport(paramCoord) ? this->EvaluateOnSupport(paramCoord) : 0.0;
}

double BasisFunction::EvaluateDerivative(int derivative, double param_coord) const {
  return derivative == 0 ? Evaluate(param_coord) :
      IsCoordinateInSupport(param_coord) ? this->EvaluateDerivativeOnSupport(derivative, param_coord) : 0.0;
}

BasisFunction::BasisFunction(const KnotVector &knot_vector, int degree, uint64_t start)
    : knotVector_(knot_vector), degree_(degree), start_of_support_(start) {}

double BasisFunction::GetKnot(uint64_t knot_position) const {
  return knotVector_.knot(knot_position);
}

uint64_t BasisFunction::GetStartOfSupport() const {
  return start_of_support_;
}

int BasisFunction::GetDegree() const {
  return degree_;
}

bool BasisFunction::IsCoordinateInSupport(double param_coord) const {
  return knotVector_.IsInKnotVectorRange(param_coord)
      && (IsCoordinateInSupportSpan(param_coord) || IsCoordinateSpecialCaseWithLastKnot(param_coord));
}

bool BasisFunction::IsCoordinateInSupportSpan(double param_coord) const {
  auto span = knotVector_.GetKnotSpan(param_coord);
  return !(span < start_of_support_ || span >= start_of_support_ + degree_ + 1);
}

bool BasisFunction::IsCoordinateSpecialCaseWithLastKnot(double param_coord) const {
  return NumericSettings<double>::AreEqual(param_coord, knotVector_.GetLastKnot()) &&
      NumericSettings<double>::AreEqual(knotVector_.GetLastKnot(), knotVector_.knot(start_of_support_ + degree_ + 1));
}
