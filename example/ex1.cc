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
#include <any>

#include "knot_vector.h"
#include "control_point.h"
#include "b_spline.h"

int main() {
  std::array<Degree, 1> degree = {Degree{2}};
  std::vector<baf::ControlPoint> control_points = {
      baf::ControlPoint(std::vector<double>({4.0, 2.0})),
      baf::ControlPoint(std::vector<double>({7.0, 1.5})),
      baf::ControlPoint(std::vector<double>({11.0, 2.0})),
      baf::ControlPoint(std::vector<double>({2.5, 4.5})),
      baf::ControlPoint(std::vector<double>({3.0, 4.3})),
      baf::ControlPoint(std::vector<double>({7.0, 3.0})),
      baf::ControlPoint(std::vector<double>({8.0, 2.5})),
      baf::ControlPoint(std::vector<double>({8.0, 1.0}))
  };
  std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_ptr = {
      std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2}, ParamCoord{3},
                                                         ParamCoord{4}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}}))};
  spl::BSpline<1> b_spline(knot_vector_ptr, degree, control_points);

  b_spline.Evaluate({ParamCoord{1.0}}, {0});

  return 0;
}
