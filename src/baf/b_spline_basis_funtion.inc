/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

// NOLINT(build/header_guard)

inline Degree BSplineBasisFunction::GetDegree() const {
  return degree_;
}

inline ParametricCoordinate BSplineBasisFunction::GetStartKnot() const {
  return start_knot_;
}

inline ParametricCoordinate BSplineBasisFunction::GetEndKnot() const {
  return end_knot_;
}

inline bool BSplineBasisFunction::IsCoordinateInSupport(ParametricCoordinate const &param_coord) const {
  return (start_knot_ <= param_coord && param_coord < end_knot_)
      || (end_knot_is_last_knot_ && param_coord == end_knot_);
}
