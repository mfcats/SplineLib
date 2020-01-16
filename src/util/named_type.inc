/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

// NOLINT(build/header_guard)


template<typename TYPE, typename NAME>
constexpr NamedType<TYPE, NAME>::NamedType(TYPE const &value) : value_(value) {}

template<typename TYPE, typename NAME>
constexpr NamedType<TYPE, NAME>::NamedType(TYPE &&value) noexcept : value_(std::move(value)) {}

template<typename TYPE, typename NAME>
constexpr NamedType<TYPE, NAME>::NamedType(NamedType<TYPE, NAME> &&other) noexcept : value_(std::move(other.value_)) {}

template<typename TYPE, typename NAME>
constexpr NamedType<TYPE, NAME> & NamedType<TYPE, NAME>::operator=(NamedType<TYPE, NAME> &&rhs) noexcept {
  value_ = std::move(rhs.value_);
  return (*this);
}

template<typename TYPE, typename NAME>
constexpr NamedType<TYPE, NAME> & NamedType<TYPE, NAME>::operator++() {
  ++value_;
  return (*this);
}

template<typename TYPE, typename NAME>
constexpr const NamedType<TYPE, NAME> NamedType<TYPE, NAME>::operator++(int) {
  NamedType<TYPE, NAME> result(*this);
  ++(*this);
  return result;
}

template<typename TYPE, typename NAME>
constexpr NamedType<TYPE, NAME> & NamedType<TYPE, NAME>::operator--() {
  --value_;
  return (*this);
}

template<typename TYPE, typename NAME>
constexpr const NamedType<TYPE, NAME> NamedType<TYPE, NAME>::operator--(int) {
  NamedType<TYPE, NAME> result(*this);
  --(*this);
  return result;
}

template<typename TYPE, typename NAME>
TYPE & NamedType<TYPE, NAME>::Get() {
  return value_;
}

template<typename TYPE, typename NAME>
constexpr TYPE const & NamedType<TYPE, NAME>::Get() const {
  return value_;
}

template<typename TYPE, typename NAME>
constexpr bool operator==(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs) {
  return numeric_settings::AreEqual<TYPE>(lhs.value_, rhs.value_);
}

template<typename TYPE, typename NAME>
constexpr bool operator!=(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs) {
  return !numeric_settings::AreEqual<TYPE>(lhs.value_, rhs.value_);
}

template<typename TYPE, typename NAME>
constexpr bool operator>(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs) {
  return (lhs.value_ > rhs.value_);
}

template<typename TYPE, typename NAME>
constexpr bool operator<(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs) {
  return (lhs.value_ < rhs.value_);
}

template<typename TYPE, typename NAME>
constexpr bool operator>=(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs) {
  if (lhs == rhs) return true;
  return (lhs > rhs);
}

template<typename TYPE, typename NAME>
constexpr bool operator<=(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs) {
  if (lhs == rhs) return true;
  return (lhs < rhs);
}

template<typename TYPE, typename NAME>
constexpr NamedType<TYPE, NAME> operator+(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs) {
  return NamedType<TYPE, NAME>{lhs.value_ + rhs.value_};
}

template<typename TYPE, typename NAME>
constexpr NamedType<TYPE, NAME> operator-(NamedType<TYPE, NAME> const &lhs, NamedType<TYPE, NAME> const &rhs) {
  return NamedType<TYPE, NAME>{lhs.value_ - rhs.value_};
}

template<typename TYPE, typename NAME>
constexpr NamedType<TYPE, NAME> operator*(TYPE factor_lhs, NamedType<TYPE, NAME> const &rhs) {
return NamedType<TYPE, NAME>{factor_lhs * rhs.value_};
}

template<typename TYPE, typename NAME>
constexpr NamedType<TYPE, NAME> operator*(NamedType<TYPE, NAME> const &lhs, TYPE factor_rhs) {
  return NamedType<TYPE, NAME>{lhs.value_ * factor_rhs};
}
