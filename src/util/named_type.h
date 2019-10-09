/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_UTIL_NAMED_TYPE_H_
#define SRC_UTIL_NAMED_TYPE_H_

#include <utility>

#include "src/util/numeric_settings.h"

namespace splinelib::src::util {
template<typename TYPE, typename NAME>
class NamedType;

template<typename TYPE, typename NAME>
constexpr bool operator==(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs);
template<typename TYPE, typename NAME>
constexpr bool operator>(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs);
template<typename TYPE, typename NAME>
constexpr bool operator<(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs);
template<typename TYPE, typename NAME>
constexpr bool operator>=(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs);
template<typename TYPE, typename NAME>
constexpr bool operator<=(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs);
template<typename TYPE, typename NAME>
constexpr NamedType<TYPE, NAME> operator+(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs);
template<typename TYPE, typename NAME>
constexpr NamedType<TYPE, NAME> operator-(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs);

// NamedType<TYPE, NAME>s are used to check for semantic of arguments (variables) at compile time. Actual definitions of
// named types are placed at the end of this file.
// Example (parametric coordinate equal to 0.5):
//   using ParametricCoordinate = splinelib::src::util::NamedType<double, struct ParametricCoordinateName>;
//   ParametricCoordinate parametric_coordinate{0.5};
//   curve.evaluate(parametric_coordinate.Get());  // evaluates 1D spline at 0.5;
template<typename TYPE, typename NAME>
class NamedType {
 public:
  NamedType() = default;
  explicit constexpr NamedType(TYPE const &value);
  explicit constexpr NamedType(TYPE &&value) noexcept;
  constexpr NamedType(NamedType const &other) = default;
  constexpr NamedType(NamedType &&other) noexcept;
  constexpr NamedType & operator=(NamedType const &rhs) = default;
  constexpr NamedType & operator=(NamedType &&rhs) noexcept;
  ~NamedType() = default;

  constexpr NamedType & operator++();
  constexpr NamedType operator++(int);
  constexpr NamedType & operator--();
  constexpr NamedType operator--(int);

  friend constexpr bool operator== <TYPE, NAME>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr bool operator> <TYPE, NAME>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr bool operator< <TYPE, NAME>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr bool operator>= <TYPE, NAME>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr bool operator<= <TYPE, NAME>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr NamedType operator+ <TYPE, NAME>(NamedType const &lhs, NamedType const &rhs);
  friend constexpr NamedType operator- <TYPE, NAME>(NamedType const &lhs, NamedType const &rhs);

  TYPE & Get();
  constexpr TYPE const & Get() const;

 private:
  TYPE value_;
};

#include "src/util/named_type.inc"
}  // namespace splinelib::src::util

namespace splinelib::src {
using Degree = util::NamedType<int, struct DegreeName>;
using Derivative = util::NamedType<int, struct DerivativeName>;
using Dimension = util::NamedType<int, struct DimensionName>;
using KnotSpan = util::NamedType<int, struct KnotSpanName>;
using Multiplicity = util::NamedType<int, struct MultiplicityName>;
using ParametricCoordinate = util::NamedType<double, struct ParametricCoordinateName>;
using Tolerance = util::NamedType<double, struct ToleranceName>;
using Weight = util::NamedType<double, struct WeightName>;
}  // namespace splinelib::src

#endif  // SRC_UTIL_NAMED_TYPE_H_
