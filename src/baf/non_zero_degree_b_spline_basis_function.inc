/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

// NOLINT(build/header_guard)

inline double NonZeroDegreeBSplineBasisFunction::ComputeLeftQuotient(ParametricCoordinate const &parametric_coordinate)
const {
  return ((parametric_coordinate - BSplineBasisFunction::GetStartKnot()).Get() * left_denominator_inverse_);
}

inline double NonZeroDegreeBSplineBasisFunction::ComputeRightQuotient(ParametricCoordinate const &parametric_coordinate)
const {
  return ((BSplineBasisFunction::GetEndKnot() - parametric_coordinate).Get() * right_denominator_inverse_);
}
