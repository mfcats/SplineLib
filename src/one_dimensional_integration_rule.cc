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
