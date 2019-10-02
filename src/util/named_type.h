/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_UTIL_NAMED_TYPE_H_
#define SRC_UTIL_NAMED_TYPE_H_

namespace splinelib::src::util {
template<typename T, typename Parameter>
class NamedType {
 public:
  NamedType() = default;

  explicit NamedType(T const &value) : value_(value) {}
  explicit NamedType(T &&value) noexcept : value_(std::move(value)) {}

  constexpr T &get() { return value_; }
  constexpr T const &get() const { return value_; }

  constexpr NamedType<T, Parameter> operator+(const NamedType<T, Parameter> &rhs) const {
    return NamedType<T, Parameter>{value_ + rhs.get()};
  }

  constexpr NamedType<T, Parameter> operator-(const NamedType<T, Parameter> &rhs) const {
    return NamedType<T, Parameter>{value_ - rhs.get()};
  }

  constexpr NamedType<T, Parameter> operator*(const NamedType<T, Parameter> &rhs) const {
    return NamedType<T, Parameter>{value_ * rhs.get()};
  }

  constexpr bool operator==(const NamedType<T, Parameter> &rhs) const {
    return value_ == rhs.get();
  }

  constexpr bool operator>(const NamedType<T, Parameter> &rhs) const {
    return value_ > rhs.get();
  }

  constexpr bool operator<(const NamedType<T, Parameter> &rhs) const {
    return value_ < rhs.get();
  }

  constexpr bool operator>=(const NamedType<T, Parameter> &rhs) const {
    return value_ >= rhs.get();
  }

  constexpr bool operator<=(const NamedType<T, Parameter> &rhs) const {
    return value_ <= rhs.get();
  }

 private:
  T value_;
};
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
