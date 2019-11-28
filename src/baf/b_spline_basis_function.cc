/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "src/baf/b_spline_basis_function.h"

#include <string>

#include "src/baf/non_zero_degree_b_spline_basis_function.h"
#include "src/baf/zero_degree_b_spline_basis_function.h"

namespace splinelib::src::baf {
BSplineBasisFunction * BSplineBasisFunction::CreateDynamic(KnotVector const &knot_vector,
                                                           KnotSpan const &start_of_support, Degree const &degree) {
  if (degree < Degree{0}) {
    throw std::logic_error("splinelib::src::baf::BSplineBasisFunction::CreateDynamic: Basis function degree must be "
                             "positive. Given degree is " + std::to_string(degree.Get()));
  }
  if (degree == Degree{0}) {
    return new ZeroDegreeBSplineBasisFunction(knot_vector, start_of_support);
  }
  return new NonZeroDegreeBSplineBasisFunction(knot_vector, start_of_support, degree);
}

double BSplineBasisFunction::Evaluate(ParametricCoordinate const &parametric_coordinate) const {
  if (IsCoordinateInSupport(parametric_coordinate)) return this->EvaluateOnSupport(parametric_coordinate);
  return 0.0;
}

double BSplineBasisFunction::EvaluateDerivative(ParametricCoordinate const &parametric_coordinate,
                                                Derivative const &derivative) const {
  if (derivative.Get() == 0) return Evaluate(parametric_coordinate);
  if (IsCoordinateInSupport(parametric_coordinate))
    return EvaluateDerivativeOnSupport(parametric_coordinate, derivative);
  return 0.0;
}

BSplineBasisFunction::BSplineBasisFunction(KnotVector const &knot_vector, KnotSpan const &start_of_support,
                                           Degree const &degree) : degree_(degree) {
  auto const start_index = start_of_support.Get();
  auto const degree_index = degree.Get();
  start_knot_ = knot_vector[start_index];
  end_knot_ = knot_vector[start_index + degree_index + 1];
  end_knot_is_last_knot_ = knot_vector.IsLastKnot(end_knot_);
}
}  // namespace splinelib::src::baf
