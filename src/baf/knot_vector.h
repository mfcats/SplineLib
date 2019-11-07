/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_BAF_KNOT_VECTOR_H_
#define SRC_BAF_KNOT_VECTOR_H_

#include <array>
#include <initializer_list>
#include <memory>
#include <vector>

#include "src/util/named_type.h"
#include "src/util/stl_container_access.h"

namespace splinelib::src::baf {
class KnotVector {
 public:
  using ConstKnotIterator = std::vector<ParametricCoordinate>::const_iterator;
  using KnotIterator = std::vector<ParametricCoordinate>::iterator;

  KnotVector() = default;
  explicit KnotVector(std::vector<ParametricCoordinate> const &knots);
  KnotVector(std::initializer_list<ParametricCoordinate> const &knots) noexcept;
  KnotVector(ConstKnotIterator begin, ConstKnotIterator end);
  KnotVector(const KnotVector &other) = default;
  KnotVector(KnotVector &&other) noexcept;
  KnotVector & operator=(KnotVector const &other) = default;
  KnotVector & operator=(KnotVector &&other) noexcept;
  // TODO(Corinna, Christoph): this should not be a constructor but implemented where it is used?
  KnotVector(std::vector<ParametricCoordinate> coords, Degree degree, int nbControlPoints);
  virtual ~KnotVector() = default;

  ParametricCoordinate const &operator[](int index) const;
  friend bool operator==(KnotVector const &lhs, KnotVector const &rhs);
  friend KnotVector operator-(KnotVector const &lhs, KnotVector const &rhs);

  ConstKnotIterator begin() const;
  ConstKnotIterator end() const;
  KnotIterator begin();
  KnotIterator end();

  virtual int GetNumberOfKnots() const;
  int GetNumberOfDifferentKnots() const;

  ParametricCoordinate GetFirstKnot() const;
  ParametricCoordinate GetLastKnot() const;
  virtual ParametricCoordinate GetKnot(int index) const;

  virtual KnotSpan GetKnotSpan(ParametricCoordinate const &parametric_coordinate) const;
  virtual Multiplicity GetMultiplicity(ParametricCoordinate const &parametric_coordinate) const;
  virtual bool IsInRange(ParametricCoordinate const &parametric_coordinate) const;
  virtual bool IsLastKnot(ParametricCoordinate const &parametric_coordinate) const;

  void InsertKnot(ParametricCoordinate const &parametric_coordinate);
  void RemoveKnot(ParametricCoordinate const &parametric_coordinate);

  // Check if absolute distance between all knots is smaller than the given tolerance.
  bool AreEqual(KnotVector const &rhs, double tolerance) const;

 private:
  std::vector<ParametricCoordinate> knots_;
};

// Check if absolute distance between all knots is smaller than the epsilon defined in NumericSettings.
bool operator==(KnotVector const &lhs, KnotVector const &rhs);
KnotVector operator-(KnotVector const &lhs, KnotVector const &rhs);

#include "src/baf/knot_vector.inc"

template<int PARAMETRIC_DIMENSIONALITY>
using KnotVectors = std::array<std::shared_ptr<KnotVector>, PARAMETRIC_DIMENSIONALITY>;
}  // namespace splinelib::src::baf

#endif  // SRC_BAF_KNOT_VECTOR_H_
