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

// KnotVectors represent a sequence U = {u_0, u_1, ..., u_m} of m+1 non-decreasing real numbers, called knots or
// parametric coordinates, on which basis functions can be defined. For a knot vector of degree p there is the
// additional condition that the first p+1 knots and the last p+1 knots are respectively equal a = u_0 = ... = u_p and
// b = u_{m-p} = ... = u_m.
// Example (knot vector of degree 2 with repeated knot 0.5 and m+1 = i knots):
//   KnotVector knot_vector({0.0, 0.0, 0.0, 0.5, 0.5, 0.75, 1.0, 1.0, 1.0});
//   sixth_knot = knot_vector[5];  // Returns sixth knot 0.75.
//   number_of_knots = knot_vector.GetNumberOfKnots();  // Returns m+1 = 9.
//   number_of_different_knots = knot_vector.GetNumberOfDifferentKnots();  // Returns 4 = |U|.
//   multiplicity_of_0_5 = knot_vector.GetMultiplicity(ParametricCoordinate{0.5});  // Returns multiplicity of u = 0.5
//   Multiplicity{2}.
//   is_1_2_in_range = knot_vector.IsInRange(ParametricCoordinate{1.2});  // Returns false as 1.2 is greater u_m = 1.0.
//   is_1_0_last_knot = knot_vector.IsLastKnot(ParametricCoordinate{1.0});  // Returns true as u_m = 1.0.
namespace splinelib::src::baf {
class KnotVector {
 public:
  using ConstKnotIterator = std::vector<ParametricCoordinate>::const_iterator;
  using KnotIterator = std::vector<ParametricCoordinate>::iterator;

  // TODO(Corinna, Christoph): remove default constructor? not possible yet because needed for mocking in zero-degree-
  //  and b-spline-basis-function-test.
  KnotVector() = default;
  explicit KnotVector(std::vector<ParametricCoordinate> knots);
  KnotVector(std::initializer_list<ParametricCoordinate> const & knots);
  KnotVector(ConstKnotIterator begin, ConstKnotIterator end);
  KnotVector(KnotVector const &other) = default;
  KnotVector(KnotVector &&other) noexcept;
  KnotVector & operator=(KnotVector const &other) = default;
  KnotVector & operator=(KnotVector &&other) noexcept;
  virtual ~KnotVector() = default;

  virtual ParametricCoordinate operator[](int index) const;
  // Check if absolute distance between all knots is smaller than the epsilon defined in NumericSettings.
  friend bool operator==(KnotVector const &lhs, KnotVector const &rhs);

  ConstKnotIterator begin() const;
  ConstKnotIterator end() const;

  virtual int GetNumberOfKnots() const;
  int GetNumberOfDifferentKnots() const;

  ParametricCoordinate GetFirstKnot() const;
  ParametricCoordinate GetLastKnot() const;
  virtual ParametricCoordinate GetKnot(int index) const;

  // The i-th knot span is defined as the half-open interval (u_i, u_{i+1}]. The last knot is defined to be in the last
  // non-zero knot span. For U = {0.0, 0.0, 0.0, 0.5, 0.5, 0.75, 1.0, 1.0, 1.0} u = 0.0 is in knot span 2, u = 0.5 in
  // knot span 4 and u = 0.75 and u = 1.0 are in knot span 5.
  // vector U = {0.0, 0.0, 0.0, 0.5, 0.5, 0.75, 1.0, 1.0, 1.0};
  virtual KnotSpan GetKnotSpan(ParametricCoordinate const &parametric_coordinate) const;
  virtual Multiplicity GetMultiplicity(ParametricCoordinate const &parametric_coordinate) const;

  // Checks if given parametric coordinate is greater equal than first knot and less equal than last knot.
  virtual bool IsInRange(ParametricCoordinate const &parametric_coordinate) const;
  virtual bool IsLastKnot(ParametricCoordinate const &parametric_coordinate) const;

  void InsertKnot(ParametricCoordinate const &parametric_coordinate);
  bool RemoveKnot(ParametricCoordinate const &parametric_coordinate);

  // Check if absolute distance between all knots is smaller than the given tolerance.
  bool AreEqual(KnotVector const &rhs, Tolerance const &tolerance) const;

 private:
  void ThrowIfKnotVectorContainsLessThanTwoKnotsOrIsNotNonDecreasing() const;
  std::vector<ParametricCoordinate> knots_;
};

// Check if absolute distance between all knots is smaller than the epsilon defined in NumericSettings.
bool operator==(KnotVector const &lhs, KnotVector const &rhs);

#include "src/baf/knot_vector.inc"

template<int PARAMETRIC_DIMENSIONALITY>
using KnotVectors = std::array<std::shared_ptr<KnotVector>, PARAMETRIC_DIMENSIONALITY>;
}  // namespace splinelib::src::baf

#endif  // SRC_BAF_KNOT_VECTOR_H_
