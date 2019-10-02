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

#include <array>
#include <initializer_list>
#include <utility>
#include <vector>

#include "named_type.h"

namespace splinelib::src::baf {
class KnotVector {
 public:
  using ConstKnotIterator = std::vector<ParametricCoordinate>::const_iterator;
  using KnotIterator = std::vector<ParametricCoordinate>::iterator;

  KnotVector() = default;
  KnotVector(const KnotVector &knotVector);
  KnotVector(KnotVector &&knotVector) noexcept;
  explicit KnotVector(std::vector<ParametricCoordinate> knots);
  KnotVector(std::initializer_list<ParametricCoordinate> knots) noexcept;
  KnotVector(ConstKnotIterator begin, ConstKnotIterator end);
  KnotVector(std::vector<ParametricCoordinate> coords, Degree degree, int nbControlPoints);

  virtual ~KnotVector() = default;

  KnotVector operator-(const KnotVector &rhs) const;
  KnotVector &operator=(const KnotVector &other);
  KnotVector &operator=(KnotVector &&other) noexcept;
  // Check if absolute distance between all knots is smaller than the epsilon defined in
  // NumericSettings.
  bool operator==(const KnotVector &rhs) const;
  bool AreEqual(const KnotVector &rhs, double tolerance) const;
  ParametricCoordinate &operator[](size_t index);

  virtual ParametricCoordinate GetKnot(size_t index) const;
  ParametricCoordinate GetLastKnot() const;
  virtual KnotSpan GetKnotSpan(ParametricCoordinate param_coord) const;
  virtual size_t GetMultiplicity(ParametricCoordinate param_coord) const;
  virtual size_t GetNumberOfKnots() const;
  int GetNumberOfDifferentKnots() const;

  ConstKnotIterator begin() const;
  ConstKnotIterator end() const;

  KnotIterator begin();
  KnotIterator end();

  virtual bool IsInKnotVectorRange(const ParametricCoordinate &param_coord) const;
  virtual bool IsLastKnot(const ParametricCoordinate &param_coord) const;

  void InsertKnot(const ParametricCoordinate &param_coord);
  void RemoveKnot(const ParametricCoordinate &param_coord);

 private:
  std::vector<ParametricCoordinate> knots_;
};

template<int PARAMETRIC_DIMENSIONALITY>
using KnotVectors = std::array<std::shared_ptr<KnotVector>, PARAMETRIC_DIMENSIONALITY>;
}  // namespace splinelib::src::baf

#endif  // SRC_BAF_KNOT_VECTOR_H_
