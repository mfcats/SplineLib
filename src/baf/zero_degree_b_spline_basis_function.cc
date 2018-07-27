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

#include <cmath>

#include "numeric_settings.h"

baf::ZeroDegreeBSplineBasisFunction::ZeroDegreeBSplineBasisFunction(const baf::KnotVector &knot_vector,
                                                                    uint64_t start_of_support) :
                                                                    BasisFunction(knot_vector, Degree{0}, KnotSpan{start_of_support}),
                                                                    start_knot_(knot_vector.GetKnot(start_of_support)),
                                                                    end_knot_(knot_vector.GetKnot(start_of_support + 1))
                                                                    {}

double baf::ZeroDegreeBSplineBasisFunction::EvaluateOnSupport(ParamCoord /* param_coord*/) const {
  return util::NumericSettings<double>::AreEqual(start_knot_.get(), end_knot_.get()) ? 0.0 : 1.0;
}

double baf::ZeroDegreeBSplineBasisFunction::EvaluateDerivativeOnSupport(ParamCoord /*param_coord*/,
                                                                        Derivative /*degree*/) const {
  return 0.0;
}
