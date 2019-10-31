/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <array>
#include <memory>
#include <vector>

#include "src/spl/b_spline.h"
#include "src/baf/control_point.h"
#include "src/baf/knot_vector.h"

using ParametricCoordinate = splinelib::src::ParametricCoordinate;

int main() {
  std::array<splinelib::src::Degree, 1> degree = {splinelib::src::Degree{2}};
  std::vector<splinelib::src::baf::ControlPoint> control_points = {
      splinelib::src::baf::ControlPoint(std::vector<double>({4.0, 2.0})),
      splinelib::src::baf::ControlPoint(std::vector<double>({7.0, 1.5})),
      splinelib::src::baf::ControlPoint(std::vector<double>({11.0, 2.0})),
      splinelib::src::baf::ControlPoint(std::vector<double>({2.5, 4.5})),
      splinelib::src::baf::ControlPoint(std::vector<double>({3.0, 4.3})),
      splinelib::src::baf::ControlPoint(std::vector<double>({7.0, 3.0})),
      splinelib::src::baf::ControlPoint(std::vector<double>({8.0, 2.5})),
      splinelib::src::baf::ControlPoint(std::vector<double>({8.0, 1.0}))
  };
  splinelib::src::baf::KnotVectors<1> knot_vector_ptr = {std::make_shared<splinelib::src::baf::KnotVector>(
      splinelib::src::baf::KnotVector({ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0},
                                       ParametricCoordinate{1},
                                       ParametricCoordinate{2}, ParametricCoordinate{3}, ParametricCoordinate{4},
                                       ParametricCoordinate{4},
                                       ParametricCoordinate{5}, ParametricCoordinate{5}, ParametricCoordinate{5}}))};
  splinelib::src::spl::BSpline<1> b_spline(knot_vector_ptr, degree, control_points);

  b_spline.Evaluate({ParametricCoordinate{1.0}}, {0});

  return 0;
}
