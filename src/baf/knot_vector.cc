/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "src/baf/knot_vector.h"

#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>

#include "src/util/numeric_settings.h"

namespace splinelib::src::baf {
KnotVector::KnotVector(std::vector<ParametricCoordinate> knots) : knots_(std::move(knots)) {}

KnotVector::KnotVector(const KnotVector &knotVector) = default;

KnotVector::KnotVector(KnotVector &&knotVector) noexcept : knots_(std::move(knotVector.knots_)) {}

KnotVector::KnotVector(std::initializer_list<ParametricCoordinate> knots) noexcept : knots_(knots) {}

KnotVector::KnotVector(ConstKnotIterator begin, ConstKnotIterator end) : knots_(
    std::vector<ParametricCoordinate>(begin, end)) {}

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

KnotVector KnotVector::operator-(const KnotVector &rhs) const {
  std::vector<ParametricCoordinate> differences(GetNumberOfKnots(), ParametricCoordinate{0.0});
  std::transform(begin(), end(), rhs.begin(), differences.begin(), std::minus<>());
  return KnotVector(differences);
}

KnotVector &KnotVector::operator=(const KnotVector &other) = default;

KnotVector &KnotVector::operator=(KnotVector &&other) noexcept {
  knots_ = std::move(other.knots_);
  return *this;
}

bool KnotVector::operator==(const KnotVector &rhs) const {
  return std::equal(this->begin(), this->end(), rhs.begin(), rhs.end(),
                    [&](ParametricCoordinate knot_a, ParametricCoordinate knot_b) {
                      return util::numeric_settings::AreEqual<double>(knot_a.Get(), knot_b.Get());
                    });
}

bool KnotVector::AreEqual(const KnotVector &rhs,
                          double tolerance = util::numeric_settings::GetEpsilon<double>()) const {
  return std::equal(this->begin(), this->end(), rhs.begin(), rhs.end(),
                    [&](ParametricCoordinate knot_a, ParametricCoordinate knot_b) {
                      return util::numeric_settings::AreEqual<double>(knot_a.Get(), knot_b.Get(), tolerance);
                    });
}

ParametricCoordinate &KnotVector::operator[](size_t index) {
#ifdef DEBUG
  return knots_.at(index);
#else
  return knots_[index];
#endif
}

ParametricCoordinate KnotVector::GetKnot(size_t index) const {
#ifdef DEBUG
  return knots_.at(index);
#else
  return knots_[index];
#endif
}

ParametricCoordinate KnotVector::GetLastKnot() const {
  return knots_.back();
}

KnotSpan KnotVector::GetKnotSpan(ParametricCoordinate param_coord) const {
  if (IsLastKnot(param_coord)) {
    return KnotSpan{static_cast<int>(std::lower_bound(knots_.begin(), knots_.end(), param_coord) - knots_.begin() - 1)};
  }
  return KnotSpan{static_cast<int>(std::upper_bound(knots_.begin(), knots_.end(), param_coord) - knots_.begin() - 1)};
}

size_t KnotVector::GetMultiplicity(ParametricCoordinate param_coord) const {
  return static_cast<size_t>(std::count(knots_.begin(), knots_.end(), param_coord));
}

size_t KnotVector::GetNumberOfKnots() const {
  return knots_.size();
}

int KnotVector::GetNumberOfDifferentKnots() const {
  int number = 1;
  for (size_t i = 1; i < knots_.size(); ++i) {
    if (knots_[i].Get() > knots_[i - 1].Get()) ++number;
  }
  return number;
}

KnotVector::ConstKnotIterator KnotVector::begin() const {
  return knots_.begin();
}

KnotVector::ConstKnotIterator KnotVector::end() const {
  return knots_.end();
}

KnotVector::KnotIterator KnotVector::begin() {
  return knots_.begin();
}

KnotVector::KnotIterator KnotVector::end() {
  return knots_.end();
}

bool KnotVector::IsInKnotVectorRange(const ParametricCoordinate &param_coord) const {
  return param_coord >= knots_.front() && param_coord <= knots_.back();
}

bool KnotVector::IsLastKnot(const ParametricCoordinate &param_coord) const {
  return util::numeric_settings::AreEqual<double>(param_coord.Get(), knots_.back().Get());
}

void KnotVector::InsertKnot(const ParametricCoordinate &param_coord) {
  KnotSpan knot_span = GetKnotSpan(param_coord);
  knots_.insert(knots_.begin() + knot_span.Get() + 1, param_coord);
}

void KnotVector::RemoveKnot(const ParametricCoordinate &param_coord) {
  KnotSpan knot_span = GetKnotSpan(param_coord);
  if (!IsLastKnot(param_coord)) {
    knots_.erase(knots_.begin() + knot_span.Get());
  } else {
    knots_.erase(knots_.end() - 1);
  }
}
}  // namespace splinelib::src::baf
