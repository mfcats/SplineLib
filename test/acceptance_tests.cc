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
#include <chrono>
#include <memory>
#include <iostream>

#include "b_spline.h"

std::unique_ptr<spl::BSpline<1>> GenerateSpline() {
  std::array<baf::KnotVector, 1> knot_vector =
      {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2}, ParamCoord{3},
                        ParamCoord{4}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}})};
  std::array<Degree, 1> degree = {Degree{2}};
  std::vector<baf::ControlPoint> control_points = {
      baf::ControlPoint(std::vector<double>({0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, 1.0})),
      baf::ControlPoint(std::vector<double>({1.0, 1.0})),
      baf::ControlPoint(std::vector<double>({1.5, 1.5})),
      baf::ControlPoint(std::vector<double>({2.0, 1.3})),
      baf::ControlPoint(std::vector<double>({3.0, 2.0})),
      baf::ControlPoint(std::vector<double>({4.0, 1.5})),
      baf::ControlPoint(std::vector<double>({4.0, 0.0}))
  };
  std::shared_ptr<std::array<baf::KnotVector, 1>> knot_vector_ptr =
      std::make_shared<std::array<baf::KnotVector, 1>>(knot_vector);
  return std::make_unique<spl::BSpline<1>>(knot_vector_ptr, degree, control_points);
}

int main() {
  auto b_spline = GenerateSpline();

  int repetitions = 20000000;
  std::chrono::system_clock::time_point before = std::chrono::system_clock::now();
  for (int i = 0; i < repetitions; i++) {
    b_spline->Evaluate({ParamCoord{i * 5.0 / repetitions}}, {0});
  }
  std::chrono::system_clock::time_point after = std::chrono::system_clock::now();
  double duration = std::chrono::duration_cast<std::chrono::milliseconds>(after - before).count();
  std::cout << "-------------------------------------------------------------------------------" << std::endl
            << "test result: Evaluating the spline " << repetitions << " times lasted " << duration << " milliseconds."
            << std::endl << "-------------------------------------------------------------------------------"
            << std::endl;
}
