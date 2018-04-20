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

KnotVector::KnotVector(const std::vector<double> &knots) : knots_(knots) {}

KnotVector::KnotVector(std::initializer_list<double> knots) : knots_(knots) {}

KnotVector::KnotVector(ConstKnotIterator begin, ConstKnotIterator end) : knots_(std::vector<double>(begin, end)) {}

bool KnotVector::operator==(const KnotVector &rhs) const {
  if (this->Size()!=rhs.Size()) return false;
  auto difference = this->knots_;
  std::transform(this->begin(), this->end(), rhs.begin(), difference.begin(), std::minus<double>());
  return !std::any_of(difference.begin(),
                      difference.end(),
                      [](double knt) { return std::fabs(knt) > NumericSettings<double>::kEpsilon(); });
}

double &KnotVector::operator[](uint64_t index) {
#ifdef DEBUG
  return knots_.at(index);
#else
  return knots_[index];
#endif
}

double KnotVector::knot(uint64_t index) const {
#ifdef DEBUG
  return knots_.at(index);
#else
  return knots_[index];
#endif
}

double KnotVector::GetLastKnot() const {
  return knots_.back();
}

int64_t KnotVector::GetKnotSpan(double param_coord) const {
  return NumericSettings<double>::AreEqual(param_coord, knots_.back()) ?
      std::lower_bound(knots_.begin(), knots_.end(), param_coord) - knots_.begin() - 1 :
      std::upper_bound(knots_.begin(), knots_.end(), param_coord) - knots_.begin() - 1;
}

ConstKnotIterator KnotVector::begin() const {
  return knots_.begin();
}

ConstKnotIterator KnotVector::end() const {
  return knots_.end();
}

bool KnotVector::IsInKnotVectorRange(double param_coord) const {
  return param_coord >= knots_.front() && param_coord <= knots_.back();
}

bool KnotVector::IsLastKnot(double param_coord) const {
  return NumericSettings<double>::AreEqual(param_coord, knots_.back());
}

uint64_t KnotVector::Size() const {
  return knots_.size();
}
