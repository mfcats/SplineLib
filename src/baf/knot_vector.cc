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
KnotVector::KnotVector(std::vector<ParametricCoordinate> knots) : knots_(std::move(knots)) {}

KnotVector::KnotVector(std::initializer_list<ParametricCoordinate> const &knots) noexcept : knots_(knots) {}

KnotVector::KnotVector(ConstKnotIterator begin, ConstKnotIterator end)
    : knots_(std::vector<ParametricCoordinate>(begin, end)) {}

KnotVector::KnotVector(KnotVector &&other) noexcept : knots_(std::move(other.knots_)) {}

KnotVector &KnotVector::operator=(KnotVector &&other) noexcept {
  knots_ = std::move(other.knots_);
  return (*this);
}

KnotVector::KnotVector(std::vector<ParametricCoordinate> coords, Degree degree, int nbControlPoints) {
  for (int i = 0; i <= degree.Get(); ++i) {
    knots_.emplace_back(ParametricCoordinate{0.0});
  }
  double curParametricCoordinate;
  for (int i = 1; i <= nbControlPoints - 1 - degree.Get(); ++i) {
    curParametricCoordinate = 0;
    for (int j = i; j < i + degree.Get(); ++j) {
      curParametricCoordinate += coords[j].Get();
    }
    curParametricCoordinate /= degree.Get();
    knots_.emplace_back(ParametricCoordinate{curParametricCoordinate});
  }
  for (int i = 0; i <= degree.Get(); ++i) {
    knots_.emplace_back(coords[coords.size() - 1]);
  }
}

int KnotVector::GetNumberOfDifferentKnots() const {
  int number_of_knot_spans = 0;
  for (int i = 1; i < static_cast<int>(knots_.size()); ++i) {
    if (knots_[i] > knots_[i - 1]) ++number_of_knot_spans;
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
  return util::numeric_settings::AreEqual(parametric_coordinate.Get(), knots_.back().Get());
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

KnotVector operator-(KnotVector const &lhs, KnotVector const &rhs) {
  std::vector<ParametricCoordinate> differences(lhs.GetNumberOfKnots(), ParametricCoordinate{0.0});
  std::transform(lhs.begin(), rhs.end(), rhs.begin(), differences.begin(), std::minus<>());
  return KnotVector(differences);
}
}  // namespace splinelib::src::baf
