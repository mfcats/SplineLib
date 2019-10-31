/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

// NOLINT(build/header_guard)

template<int PARAMETRIC_DIMENSIONALITY>
std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>> GetSpline(std::any const &spline) {
  try {
    return std::any_cast<std::shared_ptr<spl::BSpline<PARAMETRIC_DIMENSIONALITY>>>(spline);
  } catch (std::bad_any_cast &msg) {
    try {
      return std::any_cast<std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>>>(spline);
    } catch (std::bad_any_cast &msg) {
      throw std::runtime_error("Input has to be a shared pointer to a b-spline or nurbs of declared dimension");
    }
  }
}

template<int PARAMETRIC_DIMENSIONALITY>
bool IsRational(std::any const &spline) {
  try {
    std::any_cast<std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>>>(spline);
    return true;
  } catch (std::bad_any_cast &msg) {
    try {
      std::any_cast<std::shared_ptr<spl::BSpline<PARAMETRIC_DIMENSIONALITY>>>(spline);
      return false;
    } catch (std::bad_any_cast &msg) {
      throw std::runtime_error("Input has to be a shared pointer to a b-spline or nurbs of declared dimension");
    }
  }
}
