/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "zero_degree_b_spline_basis_function.h"

baf::ZeroDegBSplBasFnc::ZeroDegBSplBasFnc(const baf::KnotVector &knot_vector, const KnotSpan &start_of_support) :
    BasisFunction(knot_vector, Degree{0}, start_of_support) {}

double baf::ZeroDegBSplBasFnc::EvaluateOnSupport(const ParamCoord /*param_coord*/&) const {
  return 1.0;
}

double baf::ZeroDegBSplBasFnc::EvaluateDerivativeOnSupport(const ParamCoord /*param_coord*/&,
                                                           const Derivative /*degree*/&) const {
  return 0.0;
}
