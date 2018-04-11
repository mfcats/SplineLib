/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <cmath>

#include "numeric_settings.h"

#include "zero_degree_b_spline_basis_function.h"

ZeroDegreeBSplineBasisFunction::ZeroDegreeBSplineBasisFunction(KnotVector knot_vector,
                                                               uint64_t start_of_support)
    : BasisFunction(knot_vector, 0, start_of_support) {}

double ZeroDegreeBSplineBasisFunction::EvaluateOnSupport(double param_coord) const {
  return std::fabs(GetKnot(GetStartOfSupport()) - GetKnot(GetStartOfSupport() + 1))
             < NumericSettings<double>::kEpsilon() ? 0.0 : 1.0;
}

double ZeroDegreeBSplineBasisFunction::EvaluateDerivativeOnSupport(
    Derivative derivative, double param_coord) const {
  return 0.0;
}
