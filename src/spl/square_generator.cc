/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "square_generator.h"

spl::SquareGenerator::SquareGenerator() : degree_(2), number_of_knots_(6) {}

spl::SquareGenerator::SquareGenerator(u_int64_t degree, u_int64_t number_of_knots) : degree_(degree),
                                                                                     number_of_knots_(number_of_knots)
                                                                                     {}


std::unique_ptr<spl::BSpline<2>> spl::SquareGenerator::CreateSquare() const {
  auto parameter_space = GenerateParameterSpace();
  auto physical_space = GeneratePhysicalSpace();
  return std::make_unique<BSpline<2>>(parameter_space, physical_space);
}

spl::ParameterSpace<2> spl::SquareGenerator::GenerateParameterSpace() const {
  std::vector<ParamCoord> knots(number_of_knots_, ParamCoord{0.0});
  for (auto knot = knots.begin() + degree_ + 1; knot != knots.end() - degree_ - 1; ++knot) {
    *knot = *(knot - 1) + ParamCoord{1.0/(number_of_knots_ - 2.0 * degree_ - 1)};
  }
  std::fill(knots.end()-degree_-1, knots.end(), ParamCoord{1.0});
  baf::KnotVector knotVector(knots);
  return spl::ParameterSpace<2>({knotVector, knotVector}, {static_cast<int>(degree_), static_cast<int>(degree_)});
}

spl::PhysicalSpace<2> spl::SquareGenerator::GeneratePhysicalSpace() const {
  u_int64_t num_cps = number_of_knots_ - degree_ - 1;
  double delta = 2.0/(num_cps - 1.0);
  std::vector<double> coordinates(num_cps, - 1.0);
  double val = -1.0;
  for (double &coordinate : coordinates) {
    coordinate = val;
    val += delta;
  }
  std::vector<baf::ControlPoint> cps(num_cps * num_cps, baf::ControlPoint({0.0, 0.0}));
  size_t cp_it = 0;
  for (u_int64_t y_it = 0; y_it < num_cps; ++y_it) {
    for (u_int64_t x_it = 0; x_it < num_cps; ++x_it) {
      cps[cp_it++] = baf::ControlPoint(std::vector({coordinates[x_it], coordinates[y_it]}));
    }
  }
  return spl::PhysicalSpace<2>(cps, {static_cast<int>(num_cps), static_cast<int>(num_cps)});
}
