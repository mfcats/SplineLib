/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_NAMEDTYPE_H
#define SPLINELIB_NAMEDTYPE_H

namespace util {
template<typename T, typename Parameter>
class NamedType {
 public:
  explicit NamedType(T const &value) : value_(value) {}
  explicit NamedType(T &&value) : value_(std::move(value)) {}
  T &get() { return value_; }
  T const &get() const { return value_; }
  NamedType<T, Parameter> operator- (const NamedType<T, Parameter> &rhs) const {
    return NamedType<T, Parameter>{value_ - rhs.get()};
  }
  bool operator> (const NamedType<T, Parameter> &rhs) const {
    return value_ > rhs.get();
  }
  bool operator< (const NamedType<T, Parameter> &rhs) const {
    return value_ < rhs.get();
  }
 private:
  T value_;
};
}
#endif // SPLINELIB_NAMEDTYPE_H
