/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

// NOLINT(build/header_guard)

inline const ParametricCoordinate &KnotVector::operator[](int index) const {
  return GetValue(knots_, index);
}

inline KnotVector::ConstKnotIterator KnotVector::begin() const {
  return knots_.cbegin();
}

inline KnotVector::ConstKnotIterator KnotVector::end() const {
  return knots_.cend();
}

inline KnotVector::KnotIterator KnotVector::begin() {
  return knots_.begin();
}

inline KnotVector::KnotIterator KnotVector::end() {
  return knots_.end();
}

inline int KnotVector::GetNumberOfKnots() const {
  return knots_.size();
}

inline ParametricCoordinate KnotVector::GetFirstKnot() const {
  return knots_.front();
}

inline ParametricCoordinate KnotVector::GetLastKnot() const {
  return knots_.back();
}

inline ParametricCoordinate KnotVector::GetKnot(int index) const {
  return GetValue(knots_, index);
}

inline bool KnotVector::IsInRange(ParametricCoordinate const &parametric_coordinate) const {
  return ((parametric_coordinate >= knots_.front()) && (parametric_coordinate <= knots_.back()));
}
