/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <any>
#include <array>
#include <memory>
#include <vector>

#include "b_spline.h"
#include "control_point.h"
#include "knot_vector.h"
#include "vtk_writer.h"

void createVTKFromControlPoints(std::vector<baf::ControlPoint> control_points, std::string fileName) {
  std::array<Degree, 1> degree = {Degree{2}};
  std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_ptr = {
      std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                         ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
  spl::BSpline<1> b_spline(knot_vector_ptr, degree, control_points);
  std::shared_ptr<spl::BSpline<1>> b_spline_ptr = std::make_shared<spl::BSpline<1>>(b_spline);

  std::any b_spline_any = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_ptr);
  std::vector<std::any> splines_ = {b_spline_any};
  std::unique_ptr<io::VTKWriter> vtk_writer_;
  std::vector<std::vector<int>> scattering = {{100}};
  vtk_writer_->WriteFile(splines_, fileName, scattering);
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    throw std::runtime_error("Exactly 10 coordinates of control points are required for the splinecreation");
  }
  std::vector<baf::ControlPoint> control_points;
  control_points.emplace_back(baf::ControlPoint({-1.0, 0.0}));
  double x, y;
  for (int i = 2; i < argc; i += 2) {
    x = atof(argv[i-1]);
    y = atof(argv[i]);
    std::cout << x << ", " << y << std::endl;
    control_points.emplace_back(baf::ControlPoint({x, y}));
  }
  control_points.emplace_back(baf::ControlPoint({1.0, 0.0}));
  std::string fileName = "optimized_Form.vtk";
  createVTKFromControlPoints(control_points, fileName);
  return 0;
}
