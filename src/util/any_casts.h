/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_UTIL_ANY_CASTS_H_
#define SRC_UTIL_ANY_CASTS_H_

#include <any>
#include <memory>
#include <stdexcept>

#include "src/spl/b_spline.h"
#include "src/spl/nurbs.h"

namespace splinelib::src::util::any_casts {
template<int PARAMETRIC_DIMENSIONALITY>
std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>> GetSpline(std::any const &spline);

template<int MAXIMUM_EXPECTED_PARAMETRIC_DIMENSIONALITY>
int GetSplineDimension(std::any const &spline);

template<int PARAMETRIC_DIMENSIONALITY>
bool IsRational(std::any const &spline);

#include "src/util/any_casts.inc"
}  // namespace splinelib::src::util::any_casts

#endif  // SRC_UTIL_ANY_CASTS_H_
