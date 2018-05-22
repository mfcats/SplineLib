/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_KNOT_VECTOR_H_
#define SRC_KNOT_VECTOR_H_

#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace baf {
typedef std::vector<double>::const_iterator ConstKnotIterator;

class KnotVector {
 public:
  KnotVector() = default;

  explicit KnotVector(const std::vector<double> &knots);

  KnotVector(std::initializer_list<double> knots);

  KnotVector(ConstKnotIterator begin, ConstKnotIterator end);

  // Check if absolute distance between all knots is smaller than the epsilon defined in
  // NumericSettings.
  bool operator==(const KnotVector &rhs) const;

  double &operator[](uint64_t index);

  double knot(uint64_t index) const;

  double GetLastKnot() const;

  int64_t GetKnotSpan(double param_coord) const;

  ConstKnotIterator begin() const;

  ConstKnotIterator end() const;

  bool IsInKnotVectorRange(double param_coord) const;

  bool IsLastKnot(double param_coord) const;

  uint64_t NumberOfKnots() const;

 private:
  std::vector<double> knots_;
};
}

#endif  // SRC_KNOT_VECTOR_H_
