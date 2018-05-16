/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_NURBS_H
#define SPLINELIB_NURBS_H

#include <array>
#include <vector>

#include "spline.h"

namespace spl {
template<int DIM>
class NURBS : public Spline<DIM> {
 public:
  NURBS(const std::array<baf::KnotVector, DIM> &knot_vector,
        std::array<int, DIM> degree,
        const std::vector<baf::ControlPoint> &control_points) : Spline<DIM>(knot_vector, degree, control_points) {}
};
}

#endif //SPLINELIB_NURBS_H
