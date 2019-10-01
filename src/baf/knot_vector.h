/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_BAF_KNOT_VECTOR_H_
#define SRC_BAF_KNOT_VECTOR_H_

#include <initializer_list>
#include <utility>
#include <vector>

#include "named_type.h"

using ParamCoord = util::NamedType<double, struct ParamCoordParameter>;
using KnotSpan = util::NamedType<int, struct KnotSpanParameter>;
using Degree = util::NamedType<int, struct DegreeParameter>;

namespace baf {

class KnotVector {
 public:
  using ConstKnotIterator = std::vector<ParamCoord>::const_iterator;
  using KnotIterator = std::vector<ParamCoord>::iterator;

  KnotVector() = default;
  KnotVector(const KnotVector &knotVector);
  KnotVector(KnotVector &&knotVector) noexcept;
  explicit KnotVector(std::vector<ParamCoord> knots);
  KnotVector(std::initializer_list<ParamCoord> knots) noexcept;
  KnotVector(ConstKnotIterator begin, ConstKnotIterator end);
  KnotVector(std::vector<ParamCoord> coords, Degree degree, int nbControlPoints);

  virtual ~KnotVector() = default;

  KnotVector operator-(const KnotVector &rhs) const;
  KnotVector &operator=(const KnotVector &other);
  KnotVector &operator=(KnotVector &&other) noexcept;
  // Check if absolute distance between all knots is smaller than the epsilon defined in
  // NumericSettings.
  bool operator==(const KnotVector &rhs) const;
  bool AreEqual(const KnotVector &rhs, double tolerance) const;
  ParamCoord &operator[](size_t index);

  virtual ParamCoord GetKnot(size_t index) const;
  ParamCoord GetLastKnot() const;
  virtual KnotSpan GetKnotSpan(ParamCoord param_coord) const;
  virtual size_t GetMultiplicity(ParamCoord param_coord) const;
  virtual size_t GetNumberOfKnots() const;
  int GetNumberOfDifferentKnots() const;

  ConstKnotIterator begin() const;
  ConstKnotIterator end() const;

  KnotIterator begin();
  KnotIterator end();

  virtual bool IsInKnotVectorRange(const ParamCoord &param_coord) const;
  virtual bool IsLastKnot(const ParamCoord &param_coord) const;

  void InsertKnot(const ParamCoord &param_coord);
  void RemoveKnot(const ParamCoord &param_coord);

 private:
  std::vector<ParamCoord> knots_;
};
}  // namespace baf

#endif  // SRC_BAF_KNOT_VECTOR_H_
