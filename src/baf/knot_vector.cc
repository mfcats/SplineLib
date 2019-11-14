/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "src/baf/knot_vector.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <utility>

#include "src/util/numeric_settings.h"
#include "src/util/stl_container_access.h"

namespace splinelib::src::baf {
KnotVector::KnotVector(std::vector<ParametricCoordinate> knots) : knots_(std::move(knots)) {
  ThrowIfKnotVectorIsNotNonDecreasing();
}

KnotVector::KnotVector(std::initializer_list<ParametricCoordinate> const &knots) noexcept : knots_(knots) {
  ThrowIfKnotVectorIsNotNonDecreasing();
}

KnotVector::KnotVector(ConstKnotIterator begin, ConstKnotIterator end)
    : knots_(std::vector<ParametricCoordinate>(begin, end)) {}

KnotVector::KnotVector(KnotVector &&other) noexcept : knots_(std::move(other.knots_)) {}

KnotVector &KnotVector::operator=(KnotVector &&other) noexcept {
  knots_ = std::move(other.knots_);
  return (*this);
}

int KnotVector::GetNumberOfDifferentKnots() const {
  int number_of_knot_spans = 0;
  for (int knot_index = 1; knot_index < static_cast<int>(knots_.size()); ++knot_index) {
    if (GetValue(knots_, knot_index) > GetValue(knots_, knot_index - 1)) ++number_of_knot_spans;
  }
  return (number_of_knot_spans + 1);
}

KnotSpan KnotVector::GetKnotSpan(ParametricCoordinate const &parametric_coordinate) const {
  auto knots_begin = knots_.begin(), knots_end = knots_.end();
  if (IsLastKnot(parametric_coordinate)) {
    return KnotSpan{static_cast<int>(
        std::distance(knots_begin, std::lower_bound(knots_begin, knots_end, parametric_coordinate) - 1))};
  }
  return KnotSpan{static_cast<int>(
      std::distance(knots_begin, std::upper_bound(knots_begin, knots_end, parametric_coordinate) - 1))};
}

Multiplicity KnotVector::GetMultiplicity(ParametricCoordinate const &parametric_coordinate) const {
  return Multiplicity(std::count(knots_.begin(), knots_.end(), parametric_coordinate));
}

bool KnotVector::IsLastKnot(ParametricCoordinate const &parametric_coordinate) const {
  return util::numeric_settings::AreEqual(parametric_coordinate.Get(), GetBack(knots_).Get());
}

void KnotVector::InsertKnot(ParametricCoordinate const &parametric_coordinate) {
  auto const knot_span = GetKnotSpan(parametric_coordinate).Get();
  knots_.insert(knots_.begin() + knot_span + 1, parametric_coordinate);
}

void KnotVector::RemoveKnot(ParametricCoordinate const &parametric_coordinate) {
  auto const knot_span = GetKnotSpan(parametric_coordinate).Get();
  if (!IsLastKnot(parametric_coordinate)) {
    knots_.erase(knots_.begin() + knot_span);
  } else {
    knots_.erase(knots_.end() - 1);
  }
}

void KnotVector::ThrowIfKnotVectorIsNotNonDecreasing() const {
  for (uint64_t knot_index = 1; knot_index < knots_.size(); ++knot_index) {
    if (knots_[knot_index] < knots_[knot_index - 1])
      throw std::invalid_argument("splinelib::src::baf::KnotVector::ThrowIfInvalidKnotVector: Knot vector has to be "
                                  "a sequence of non-decreasing parametric coordinates!");
  }
}

bool KnotVector::AreEqual(KnotVector const &rhs,
                          Tolerance const &tolerance = Tolerance{util::numeric_settings::GetEpsilon<double>()}) const {
  return std::equal(this->begin(), this->end(), rhs.begin(), rhs.end(),
                    [&](ParametricCoordinate knot_a, ParametricCoordinate knot_b) {
                      return util::numeric_settings::AreEqual<double>(knot_a.Get(), knot_b.Get(), tolerance.Get());
                    });
}

bool operator==(KnotVector const &lhs, KnotVector const &rhs) {
  return std::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(),
                    [&](ParametricCoordinate knot_lhs, ParametricCoordinate knot_rhs) {
                      return util::numeric_settings::AreEqual<double>(knot_lhs.Get(), knot_rhs.Get());
                    });
}
}  // namespace splinelib::src::baf
