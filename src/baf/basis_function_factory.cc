/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "basis_function_factory.h"

#include "b_spline_basis_function.h"
#include "zero_degree_b_spline_basis_function.h"

baf::BasisFunction *baf::BasisFunctionFactory::CreateDynamic(const KnotVector &knot_vector,
                                                             uint64_t start_of_support,
                                                             int degree) const {
  if (degree < 0) throw std::runtime_error("Basis function degree must be positive.");
  if (degree == 0) return new baf::ZeroDegreeBSplineBasisFunction(knot_vector, start_of_support);
  return new baf::BSplineBasisFunction(knot_vector, degree, start_of_support);
}
