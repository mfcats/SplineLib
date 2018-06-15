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

#include <cmath>
#include <algorithm>
#include <functional>

#include "numeric_settings.h"

baf::KnotVector::KnotVector(const std::vector<ParamCoord> &knots) : knots_(knots) {}

baf::KnotVector::KnotVector(const baf::KnotVector &knotVector) : knots_(knotVector.knots_) {}

baf::KnotVector::KnotVector(const baf::KnotVector &&knotVector) : knots_(std::move(knotVector.knots_)) {}

baf::KnotVector::KnotVector(std::initializer_list<ParamCoord> knots) : knots_(knots) {}

baf::KnotVector::KnotVector(ConstKnotIterator begin, ConstKnotIterator end) : knots_(std::vector<ParamCoord>(begin,
                                                                                                             end)) {}

baf::KnotVector baf::KnotVector::operator-(const baf::KnotVector &rhs) const {
  std::vector<ParamCoord> differences;
  for(int knot = 0; knot < this->GetNumberOfKnots(); knot++) {
    differences.push_back(ParamCoord{knots_[knot] - rhs.knots_[knot]});
  }
  return baf::KnotVector(differences);
}

baf::KnotVector &baf::KnotVector::operator=(const baf::KnotVector &other) {
  knots_ = other.knots_;
  return *this;
}

baf::KnotVector &baf::KnotVector::operator=(const baf::KnotVector &&other) {
  knots_ = std::move(other.knots_);
  return *this;
}

bool baf::KnotVector::operator==(const KnotVector &rhs) const {
  return std::equal(this->begin(),
                    this->end(),
                    rhs.begin(),
                    rhs.end(),
                    [&](ParamCoord a, ParamCoord b) {
                      return util::NumericSettings<double>::AreEqual(a.get(),
                                                                     b.get());
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

ParamCoord baf::KnotVector::GetLastKnot() const {
  return knots_.back();
}

u_int64_t baf::KnotVector::GetKnotSpan(ParamCoord param_coord) const {
  return static_cast<u_int64_t>(util::NumericSettings<double>::AreEqual(param_coord.get(), knots_.back().get()) ?
      std::lower_bound(knots_.begin(), knots_.end(), param_coord) - knots_.begin() - 1 :
      std::upper_bound(knots_.begin(), knots_.end(), param_coord) - knots_.begin() - 1);
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

bool baf::KnotVector::IsInKnotVectorRange(ParamCoord param_coord) const {
  return param_coord.get() >= knots_.front().get() && param_coord.get() <= knots_.back().get();
}

bool baf::KnotVector::IsLastKnot(ParamCoord param_coord) const {
  return util::NumericSettings<double>::AreEqual(param_coord.get(), knots_.back().get());
}

size_t baf::KnotVector::GetNumberOfKnots() const {
  return knots_.size();
}
