/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University
This file is part of SplineLib.
SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.
SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "src/spl/square_generator.h"

namespace splinelib::src::spl {
SquareGenerator::SquareGenerator() : degree_(Degree{2}), number_of_knots_(6) {}

SquareGenerator::SquareGenerator(Degree degree, u_int64_t number_of_knots) : degree_(std::move(degree)),
                                                                             number_of_knots_(number_of_knots) {}

std::unique_ptr<BSpline<2>> SquareGenerator::CreateSquare() const {
  auto parameter_space = std::make_shared<ParameterSpace<2>>(GenerateParameterSpace());
  auto physical_space = std::make_shared<PhysicalSpace<2>>(GeneratePhysicalSpace());
  BSplineGenerator<2> spline_generator(physical_space, parameter_space);
  return std::make_unique<BSpline<2>>(spline_generator);
}

ParameterSpace<2> SquareGenerator::GenerateParameterSpace() const {
  std::vector<ParametricCoordinate> knots(number_of_knots_, zero_);
  for (auto knot = knots.begin() + degree_.Get() + 1; knot != knots.end() - degree_.Get() - 1; ++knot) {
    *knot = *(knot - 1) + ParametricCoordinate{1.0 / (number_of_knots_ - 2.0 * degree_.Get() - 1)};
  }
  std::fill(knots.end() - degree_.Get() - 1, knots.end(), one_);
  std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vectors = {
      std::make_shared<baf::KnotVector>(knots),
      std::make_shared<baf::KnotVector>(knots)
  };
  return ParameterSpace<2>(knot_vectors, {degree_, degree_});
}

PhysicalSpace<2> SquareGenerator::GeneratePhysicalSpace() const {
  u_int64_t num_cps = number_of_knots_ - degree_.Get() - 1;
  double delta = 2.0 / (num_cps - 1.0);
  std::vector<double> coordinates(num_cps, -1.0);
  double val = -1.0;
  for (auto &coordinate : coordinates) {
    coordinate = val;
    val += delta;
  }
  std::vector<spl::ControlPoint> cps(num_cps * num_cps, spl::ControlPoint({0.0, 0.0}));
  size_t cp_it = 0;
  for (u_int64_t y_it = 0; y_it < num_cps; ++y_it) {
    for (u_int64_t x_it = 0; x_it < num_cps; ++x_it) {
      cps[cp_it++] = spl::ControlPoint(std::vector({coordinates[x_it], coordinates[y_it]}));
    }
  }
  return PhysicalSpace<2>(cps, {static_cast<int>(num_cps), static_cast<int>(num_cps)});
}
}  // namespace splinelib::src::spl
