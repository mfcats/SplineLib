/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "control_point.h"

ControlPoint::ControlPoint(std::initializer_list<double> coordinates) : coordinates_(coordinates) {}

ControlPoint::ControlPoint(const std::vector<double> &coordinates) : coordinates_(coordinates) {}

int ControlPoint::GetDimension() const {
  return static_cast<int>(coordinates_.size());
}

double ControlPoint::GetValue(int dimension) const {
#ifdef DEBUG
  return coordinates_.at(static_cast<unsigned long>(dimension));
#else
  return coordinates_[dimension];
#endif
}
