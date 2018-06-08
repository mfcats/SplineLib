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

baf::KnotVector::KnotVector(std::initializer_list<ParamCoord> knots) : knots_(knots) {}

baf::KnotVector::KnotVector(ConstKnotIterator begin, ConstKnotIterator end) : knots_(std::vector<ParamCoord>(begin, end)) {}

bool baf::KnotVector::operator==(const KnotVector &rhs) const {
  if (this->NumberOfKnots() != rhs.NumberOfKnots()) return false;
  auto difference = this->knots_;
  std::transform(this->begin(), this->end(), rhs.begin(), difference.begin(), std::minus<ParamCoord>());
  return !std::any_of(difference.begin(),
                      difference.end(),
                      [](ParamCoord knt) { return std::fabs(knt.get()) > util::NumericSettings<double>::kEpsilon(); });
}

ParamCoord &baf::KnotVector::operator[](uint64_t index) {
#ifdef DEBUG
  return knots_.at(index);
#else
  return knots_[index];
#endif
}

ParamCoord baf::KnotVector::knot(uint64_t index) const {
#ifdef DEBUG
  return knots_.at(index);
#else
  return knots_[index];
#endif
}

ParamCoord baf::KnotVector::GetLastKnot() const {
  return knots_.back();
}

int64_t baf::KnotVector::GetKnotSpan(ParamCoord param_coord) const {
  return util::NumericSettings<double>::AreEqual(param_coord.get(), knots_.back().get()) ?
         std::lower_bound(knots_.begin(), knots_.end(), param_coord) - knots_.begin() - 1 :
         std::upper_bound(knots_.begin(), knots_.end(), param_coord) - knots_.begin() - 1;
}

baf::ConstKnotIterator baf::KnotVector::begin() const {
  return knots_.begin();
}

baf::ConstKnotIterator baf::KnotVector::end() const {
  return knots_.end();
}

bool baf::KnotVector::IsInKnotVectorRange(ParamCoord param_coord) const {
  return param_coord.get() >= knots_.front().get() && param_coord.get() <= knots_.back().get();
}

bool baf::KnotVector::IsLastKnot(ParamCoord param_coord) const {
  return util::NumericSettings<double>::AreEqual(param_coord.get(), knots_.back().get());
}

uint64_t baf::KnotVector::NumberOfKnots() const {
  return knots_.size();
}
