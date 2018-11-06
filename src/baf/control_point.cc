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

baf::ControlPoint::ControlPoint(std::initializer_list<double> coordinates) : coordinates_(coordinates) {}

baf::ControlPoint::ControlPoint(std::vector<double> coordinates) : coordinates_(std::move(coordinates)) {}

int baf::ControlPoint::GetDimension() const {
  return static_cast<int>(coordinates_.size());
}

baf::ControlPoint baf::ControlPoint::operator+(const baf::ControlPoint &control_point) const {
  std::vector<double> coordinates_new;
  coordinates_new.reserve(this->GetDimension());
  for (int i = 0; i < this->GetDimension(); ++i) {
    coordinates_new.push_back(this->GetValue(i) + control_point.GetValue(i));
  }
  return ControlPoint(coordinates_new);
}

baf::ControlPoint baf::ControlPoint::Transform(const baf::ControlPoint &controlPoint,
    std::array<double, 3> x, std::array<double, 3> y, std::array<double, 3>, std::array<double, 3> o) const {
  std::vector<double> coordinates_new = {-o[0], -o[1], -o[2]};
  for (int i = 0; i < 3; ++i) {
    coordinates_new[0] += x[i] * controlPoint.GetValue(i);
    coordinates_new[1] += y[i] * controlPoint.GetValue(i);
    coordinates_new[2] += z[i] * controlPoint.GetValue(i);
  }
  return ControlPoint(coordinates_new);
}

double baf::ControlPoint::GetValue(int dimension) const {
#ifdef DEBUG
  return coordinates_.at(dimension);
#else
  return coordinates_[dimension];
#endif
}
