/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_UTIL_NUMERIC_SETTINGS_H_
#define SRC_UTIL_NUMERIC_SETTINGS_H_

#include <cmath>

namespace splinelib::src::util::numeric_settings {
template<typename TYPE>
constexpr TYPE GetEpsilon();
template<typename TYPE>
constexpr bool AreEqual(TYPE const &a, TYPE const &b, TYPE const &tolerance = GetEpsilon<TYPE>());

static constexpr double kEpsilonFactor_ = 1e1;

#include "src/util/numeric_settings.inc"
}  // namespace splinelib::src::util::numeric_settings

#endif  // SRC_UTIL_NUMERIC_SETTINGS_H_
