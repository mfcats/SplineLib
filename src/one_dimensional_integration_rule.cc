/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "one_dimensional_integration_rule.h"

OneDimensionalIntegrationRule::OneDimensionalIntegrationRule(int points)
    : number_of_points_(points) {
  switch (points) {
    case 1:points_ = {0};
      weights_ = {2};
      break;
    case 2:points_ = {-sqrt(1.0 / 3), sqrt(1.0 / 3.0)};
      weights_ = {1, 1};
      break;
    case 3:points_ = {-sqrt(3.0 / 5), 0, sqrt(3.0 / 5)};
      weights_ = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
      break;
    case 4:
      points_ = {-sqrt(3.0 / 7 + 2.0 / 7 * sqrt(6.0 / 5)),
                 -sqrt(3.0 / 7 - 2.0 / 7 * sqrt(6.0 / 5)),
                 sqrt(3.0 / 7 - 2.0 / 7 * sqrt(6.0 / 5)),
                 sqrt(3.0 / 7 + 2.0 / 7 * sqrt(6.0 / 5))};
      weights_ = {(18.0 - sqrt(30)) / 36, (18.0 + sqrt(30)) / 36,
                  (18.0 + sqrt(30)) / 36, (18.0 - sqrt(30)) / 36};
      break;
    case 5:
      points_ = {-(1.0 / 3) * sqrt(5 + 2.0 * sqrt(10.0 / 7)),
                 -(1.0 / 3) * sqrt(5 - 2.0 * sqrt(10.0 / 7)), 0,
                 (1.0 / 3) * sqrt(5 - 2.0 * sqrt(10.0 / 7)),
                 (1.0 / 3) * sqrt(5 + 2.0 * sqrt(10.0 / 7))};
      weights_ = {(322.0 - 13 * sqrt(70)) / 900, (322.0 + 13 * sqrt(70)) / 900, 128.0 / 225,
                  (322.0 + 13 * sqrt(70)) / 900, (322.0 - 13 * sqrt(70)) / 900};
      break;
    default:throw std::runtime_error("Integration rules are only implemented for up to 5 points.");
  }
}

int OneDimensionalIntegrationRule::points() const {
  return number_of_points_;
}

double OneDimensionalIntegrationRule::point(int point) const {
#ifdef DEBUG
  return points_.at(point);
#else
  return points_[point];
#endif
}

double OneDimensionalIntegrationRule::weight(int point) const {
#ifdef DEBUG
  return weights_.at(point);
#else
  return weights_[point];
#endif
}
