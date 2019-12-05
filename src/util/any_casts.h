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

#include "src/spl/b_spline.h"
#include "src/spl/nurbs.h"

// The namespace any_casts handles std::any objects that contain shared pointer to splines. If it is not known which
// parametric dimensionality a spline has or if it is rational, this class converts the std::any object correctly into a
// pointer to a spline without throwing an exception. An exception only is thrown when the std::any object does not
// contain a pointer to a spline or differs from given (maximum expected) parametric dimensionality.
// Example (std::any object containing a pointer to a two-dimensional NURBS):
//   std::any any_object = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbs_2d_pointer);
//   is_rational = IsRational<2>(any_object);  // Returns true;
//   spline_dimension = GetSplineDimension<2>(any_object);  // Returns 2.
//   spline_dimension = GetSplineDimension<1>(any_object);  // Throws std::logic_error as maximum expected dimension of
//                                                          // 1 is lower than actual dimension of 2.
//   spline_pointer = GetSpline<2>(any_object);  // Returns nurbs_2d_pointer.
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
