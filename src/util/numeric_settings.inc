/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

// NOLINT(build/header_guard)

#include <limits>

template<typename TYPE>
constexpr TYPE GetEpsilon() {
  return (kEpsilonFactor_ * std::numeric_limits<TYPE>::epsilon());
}

template<typename TYPE>
constexpr bool AreEqual(TYPE const &a, TYPE const &b, TYPE const &tolerance) {
  return (std::abs(a - b) <= tolerance);
}
