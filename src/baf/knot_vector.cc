/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "knot_vector.h"

#include <algorithm>
#include <functional>

#include "numeric_settings.h"

baf::KnotVector::KnotVector(std::vector<ParamCoord> knots) : knots_(std::move(knots)) {}

baf::KnotVector::KnotVector(const baf::KnotVector &knotVector) = default;

baf::KnotVector::KnotVector(baf::KnotVector &&knotVector) noexcept : knots_(std::move(knotVector.knots_)) {}

baf::KnotVector::KnotVector(std::initializer_list<ParamCoord> knots) noexcept : knots_(knots) {}

baf::KnotVector::KnotVector(ConstKnotIterator begin, ConstKnotIterator end) : knots_(std::vector<ParamCoord>(begin,
                                                                                                             end)) {}

baf::KnotVector baf::KnotVector::operator-(const baf::KnotVector &rhs) const {
  std::vector<ParamCoord> differences(GetNumberOfKnots(), ParamCoord{0.0});
  std::transform(begin(), end(), rhs.begin(), differences.begin(), std::minus<>());
  return baf::KnotVector(differences);
}

baf::KnotVector &baf::KnotVector::operator=(const baf::KnotVector &other) = default;

baf::KnotVector &baf::KnotVector::operator=(KnotVector &&other) noexcept {
  knots_ = std::move(other.knots_);
  return *this;
}

bool baf::KnotVector::operator==(const KnotVector &rhs) const {
  return std::equal(this->begin(),
                    this->end(),
                    rhs.begin(),
                    rhs.end(),
                    [&](ParamCoord knot_a, ParamCoord knot_b) {
                      return util::NumericSettings<double>::AreEqual(knot_a.get(), knot_b.get());
                    });
}

ParamCoord &baf::KnotVector::operator[](size_t index) {
#ifdef DEBUG
  return knots_.at(index);
#else
  return knots_[index];
#endif
}

ParamCoord baf::KnotVector::GetKnot(size_t index) const {
#ifdef DEBUG
  return knots_.at(index);
#else
  return knots_[index];
#endif
}

std::vector<ParamCoord> baf::KnotVector::GetKnots() const {
  return knots_;
}

ParamCoord baf::KnotVector::GetLastKnot() const {
  return knots_.back();
}

KnotSpan baf::KnotVector::GetKnotSpan(ParamCoord param_coord) const {
  if (IsLastKnot(param_coord)) {
    return KnotSpan{static_cast<int>(std::lower_bound(knots_.begin(), knots_.end(), param_coord) - knots_.begin() - 1)};
  }
  return KnotSpan{static_cast<int>(std::upper_bound(knots_.begin(), knots_.end(), param_coord) - knots_.begin() - 1)};
}

size_t baf::KnotVector::GetMultiplicity(ParamCoord param_coord) const {
  return static_cast<size_t>(std::count(knots_.begin(), knots_.end(), param_coord));
}

size_t baf::KnotVector::GetNumberOfKnots() const {
  return knots_.size();
}

baf::KnotVector::ConstKnotIterator baf::KnotVector::begin() const {
  return knots_.begin();
}

baf::KnotVector::ConstKnotIterator baf::KnotVector::end() const {
  return knots_.end();
}

baf::KnotVector::KnotIterator baf::KnotVector::begin() {
  return knots_.begin();
}

baf::KnotVector::KnotIterator baf::KnotVector::end() {
  return knots_.end();
}

bool baf::KnotVector::IsInKnotVectorRange(const ParamCoord &param_coord) const {
  return param_coord >= knots_.front() && param_coord <= knots_.back();
}

bool baf::KnotVector::IsLastKnot(const ParamCoord &param_coord) const {
  return util::NumericSettings<double>::AreEqual(param_coord.get(), knots_.back().get());
}

void baf::KnotVector::InsertKnot(const ParamCoord &param_coord) {
  KnotSpan knot_span = GetKnotSpan(param_coord);
  knots_.insert(knots_.begin() + knot_span.get() + 1, param_coord);
}
