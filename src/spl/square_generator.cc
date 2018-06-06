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

spl::SquareGenerator::SquareGenerator() : knot_vectors_(
    {baf::KnotVector{{ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1},
                      ParamCoord{1}}/*zero_, zero_, zero_, one_, one_, one_*/}, baf::KnotVector{
        {ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1},
         ParamCoord{1}}/*zero_, zero_, zero_, one_, one_, one_*/}}),
                                          degrees_({2, 2}),
                                          control_points_({baf::ControlPoint(std::vector<double>({-1.0, -1.0})),
                                                           baf::ControlPoint(std::vector<double>({0.0, -1.0})),
                                                           baf::ControlPoint(std::vector<double>({1.0, -1.0})),
                                                           baf::ControlPoint(std::vector<double>({-1.0, 0.0})),
                                                           baf::ControlPoint(std::vector<double>({0.0, 0.0})),
                                                           baf::ControlPoint(std::vector<double>({1.0, 0.0})),
                                                           baf::ControlPoint(std::vector<double>({-1.0, 1.0})),
                                                           baf::ControlPoint(std::vector<double>({0.0, 1.0})),
                                                           baf::ControlPoint(std::vector<double>({1.0, 1.0}))}) {}

spl::SquareGenerator::SquareGenerator(int degree, int number_of_knots) : degrees_({degree, degree}) {
  std::vector<ParamCoord > knots;
  for (int i = 0; i <= degree; i++) {
    knots.push_back(zero_);
  }
  for (double i = 1; i <= number_of_knots - 2 * degree - 2; i++) {
    knots.push_back(ParamCoord{i / (number_of_knots - 2 * degree - 1)});
  }
  for (int i = 0; i <= degree; i++) {
    knots.push_back(one_);
  }
  knot_vectors_ = {baf::KnotVector(knots), baf::KnotVector(knots)};

  for (int dimension1 = 0; dimension1 < number_of_knots - degree - 1; dimension1++) {
    for (int dimension2 = 0; dimension2 < number_of_knots - degree - 1; dimension2++) {
      control_points_.push_back(baf::ControlPoint({2.0 * dimension2 / (number_of_knots - degree - 2) - 1,
                                                   2.0 * dimension1 / (number_of_knots - degree - 2) - 1}));
    }
  }
}

std::unique_ptr<spl::BSpline<2>> spl::SquareGenerator::CreateSquare() const {
  return std::make_unique<BSpline<2>>(knot_vectors_, degrees_, control_points_);
}
