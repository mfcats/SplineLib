/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "src/baf/zero_degree_b_spline_basis_function.h"

namespace splinelib::src::baf {
ZeroDegreeBSplineBasisFunction::ZeroDegreeBSplineBasisFunction(const KnotVector &knot_vector,
                                                               const KnotSpan &start_of_support) :
    BSplineBasisFunction(knot_vector, Degree{0}, start_of_support) {}

double ZeroDegreeBSplineBasisFunction::EvaluateOnSupport(const ParametricCoordinate /*param_coord*/&) const {
  return 1.0;
}

double ZeroDegreeBSplineBasisFunction::EvaluateDerivativeOnSupport(const ParametricCoordinate /*param_coord*/&,
                                                                   const Derivative /*degree*/&) const {
  return 0.0;
}
}  // namespace splinelib::src::baf
