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

SquareGenerator::SquareGenerator() : knot_vectors_({KnotVector{0, 0, 0, 1, 1, 1}, KnotVector{0, 0, 0, 1, 1, 1}}),
                                     degrees_({2, 2}),
                                     control_points_({ControlPoint(std::vector<double>({-1.0, -1.0})),
                                                      ControlPoint(std::vector<double>({0.0, -1.0})),
                                                      ControlPoint(std::vector<double>({1.0, -1.0})),
                                                      ControlPoint(std::vector<double>({-1.0, 0.0})),
                                                      ControlPoint(std::vector<double>({0.0, 0.0})),
                                                      ControlPoint(std::vector<double>({1.0, 0.0})),
                                                      ControlPoint(std::vector<double>({-1.0, 1.0})),
                                                      ControlPoint(std::vector<double>({0.0, 1.0})),
                                                      ControlPoint(std::vector<double>({1.0, 1.0}))}) {}

SquareGenerator::SquareGenerator(int degree, int number_of_knots) : degrees_({degree, degree}) {
  std::vector<double> knots;
  for (int i = 0; i <= degree; i++) {
    knots.push_back(0);
  }
  for (double i = 1; i <= number_of_knots - 2 * degree - 2; i++) {
    knots.push_back(i / (number_of_knots - 2 * degree - 1));
  }
  for (int i = 0; i <= degree; i++) {
    knots.push_back(1);
  }
  knot_vectors_ = {KnotVector(knots), KnotVector(knots)};

  for (int dimension1 = 0; dimension1 < number_of_knots - degree - 1; dimension1++) {
    for (int dimension2 = 0; dimension2 < number_of_knots - degree - 1; dimension2++) {
      control_points_.push_back(ControlPoint({2.0 * dimension2 / (number_of_knots - degree - 2) - 1,
                                              2.0 * dimension1 / (number_of_knots - degree - 2) - 1}));
    }
  }
}

std::unique_ptr<BSpline<2>> SquareGenerator::CreateSquare() const {
  return std::make_unique<BSpline<2>>(knot_vectors_, degrees_, control_points_);
}
