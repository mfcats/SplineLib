/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_UTIL_ANY_CASTS_H_
#define SRC_UTIL_ANY_CASTS_H_

#include <any>
#include <memory>

#include "b_spline.h"
#include "nurbs.h"

namespace util {
template <int DIM>
class AnyCasts {
 public:
  static std::shared_ptr<spl::Spline<DIM>> GetSpline(std::any spline) {
    try {
      return std::any_cast<std::shared_ptr<spl::BSpline<DIM>>>(spline);
    } catch (std::bad_any_cast &msg) {
      try {
        return std::any_cast<std::shared_ptr<spl::NURBS<DIM>>>(spline);
      } catch (std::bad_any_cast &msg) {
        throw std::runtime_error("Input has to be a pointer to a b-spline or nurbs.");
      }
    }
  }

  static bool IsRational(std::any spline) {
    try {
      std::any_cast<std::shared_ptr<spl::NURBS<DIM>>>(spline);
      return true;
    } catch (std::bad_any_cast &msg) {
      return false;
    }
  }
};
}  // namespace util

#endif  // SRC_UTIL_ANY_CASTS_H_
