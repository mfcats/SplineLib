/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "src/util/any_casts.h"

namespace splinelib::src::util::any_casts {
int GetSplineDimension(const std::any &spline) {
  try {
    GetSpline<1>(spline);
    return 1;
  } catch (std::runtime_error &e) {
    try {
      GetSpline<2>(spline);
      return 2;
    } catch (std::runtime_error &e) {
      try {
        GetSpline<3>(spline);
        return 3;
      } catch (std::runtime_error &e) {
        try {
          GetSpline<4>(spline);
          return 4;
        } catch (std::runtime_error &e) {
          throw std::runtime_error(
              "Input has to be a shared pointer to a b-spline or nurbs of dimension 1, 2, 3 or 4.");
        }
      }
    }
  }
}
}  // namespace splinelib::src::util::any_casts
