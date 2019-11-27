/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

// NOLINT(build/header_guard)

inline double ControlPoint::operator[](Dimension const &dimension) const {
  return GetValue(coordinates_, dimension);
}

inline double ControlPoint::GetValueForDimension(Dimension const &dimension) const {
  return GetValue(coordinates_, dimension);
}

inline int ControlPoint::GetDimensionality() const {
  return coordinates_.size();
}

inline void ControlPoint::SetValue(Dimension dimension, double value) {
  GetValue(coordinates_, dimension) = value;
}
