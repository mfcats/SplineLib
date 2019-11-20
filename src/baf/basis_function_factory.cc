/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "src/baf/basis_function_factory.h"

#include <string>

#include "src/baf/non_zero_degree_b_spline_basis_function.h"
#include "src/baf/zero_degree_b_spline_basis_function.h"

namespace splinelib::src::baf {
BSplineBasisBasisFunction *BasisFunctionFactory::CreateDynamic(const KnotVector &knot_vector,
                                                               const KnotSpan &start_of_support,
                                                               const Degree &degree) {
  if (degree < Degree{0}) {
    std::string msg;
    msg = "Basis function degree must be positive. Given degree is " + std::to_string(degree.Get());
    throw std::runtime_error(msg);
  }
  if (degree == Degree{0}) {
    return new ZeroDegreeBSplineBasisFunction(knot_vector, start_of_support);
  }
  return new NonZeroDegreeBSplineBasisFunction(knot_vector, degree, start_of_support);
}
}  // namespace splinelib::src::baf
