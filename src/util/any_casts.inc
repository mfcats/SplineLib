/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

// NOLINT(build/header_guard)

#include <any>
#include <memory>

#include "src/spl/b_spline.h"
#include "src/spl/nurbs.h"

template<int PARAMETRIC_DIMENSIONALITY>
inline std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>> GetSpline(std::any const &spline) {
  try {
    return std::any_cast<std::shared_ptr<spl::BSpline<PARAMETRIC_DIMENSIONALITY>>>(spline);
  } catch (std::bad_any_cast &message) {
    try {
      return std::any_cast<std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>>>(spline);
    } catch (std::bad_any_cast &message) {
      throw std::invalid_argument("splinelib::src::util::any_casts::GetSpline: Input has to be a shared pointer to a "
                                  "B-spline or NURBS of dimension" + PARAMETRIC_DIMENSIONALITY);
    }
  }
}

template<int MAXIMUM_EXPECTED_PARAMETRIC_DIMENSIONALITY>
inline int GetSplineDimension(std::any const &spline) {
  try {
    GetSpline<MAXIMUM_EXPECTED_PARAMETRIC_DIMENSIONALITY>(spline);
    return MAXIMUM_EXPECTED_PARAMETRIC_DIMENSIONALITY;
  } catch (std::invalid_argument &message) {
    return GetSplineDimension<MAXIMUM_EXPECTED_PARAMETRIC_DIMENSIONALITY - 1>(spline);
  }
}
template<>
inline int GetSplineDimension<0>(std::any const &/*spline*/) {
  throw std::logic_error("splinelib::src::util::any_casts::GetSplineDimension: Parametric dimensionality smaller one"
                         " or larger than the chosen maximum encountered or the input is no shared pointer to a "
                         "b-spline or nurbs.");
}

template<int PARAMETRIC_DIMENSIONALITY>
inline bool IsRational(std::any const &spline) {
  try {
    std::any_cast<std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>>>(spline);
    return true;
  } catch (std::bad_any_cast &message) {
    try {
      std::any_cast<std::shared_ptr<spl::BSpline<PARAMETRIC_DIMENSIONALITY>>>(spline);
      return false;
    } catch (std::bad_any_cast &message) {
      throw std::invalid_argument("splinelib::src::util::any_casts::IsRational: Input has to be a shared pointer to a "
                                  "b-spline or nurbs of declared dimension");
    }
  }
}
