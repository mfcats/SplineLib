/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "src/baf/zero_degree_b_spline_basis_function.h"

#include <utility>

namespace splinelib::src::baf {
ZeroDegreeBSplineBasisFunction::ZeroDegreeBSplineBasisFunction(KnotVector const &knot_vector,
                                                               KnotSpan const &start_of_support)
    : BSplineBasisFunction(knot_vector, start_of_support, Degree{0}) {}

ZeroDegreeBSplineBasisFunction::ZeroDegreeBSplineBasisFunction(ZeroDegreeBSplineBasisFunction &&rhs) noexcept
    : BSplineBasisFunction(std::move(rhs)) {}

ZeroDegreeBSplineBasisFunction & ZeroDegreeBSplineBasisFunction::operator=(
    ZeroDegreeBSplineBasisFunction &&rhs) noexcept {
  BSplineBasisFunction::operator=(std::move(rhs));
  return (*this);
}
}  // namespace splinelib::src::baf
