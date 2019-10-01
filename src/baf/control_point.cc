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

baf::ControlPoint::ControlPoint(uint64_t dimension) : coordinates_(std::vector(dimension, 0.0)) {}

int baf::ControlPoint::GetDimension() const {
  return static_cast<int>(coordinates_.size());
}

double baf::ControlPoint::GetValue(int dimension) const {
#ifdef DEBUG
  return coordinates_.at(dimension);
#else
  return coordinates_[dimension];
#endif
}

void baf::ControlPoint::SetValue(int dimension, double value) {
  coordinates_[dimension] = value;
}

baf::ControlPoint baf::ControlPoint::operator+(const baf::ControlPoint &control_point) const {
  std::vector<double> coordinates_new;
  coordinates_new.reserve(this->GetDimension());
  for (int i = 0; i < this->GetDimension(); ++i) {
    coordinates_new.push_back(this->GetValue(i) + control_point.GetValue(i));
  }
  return ControlPoint(coordinates_new);
}

baf::ControlPoint baf::ControlPoint::operator-(const baf::ControlPoint &control_point) const {
  std::vector<double> coordinates_new;
  coordinates_new.reserve(this->GetDimension());
  for (int i = 0; i < this->GetDimension(); ++i) {
    coordinates_new.push_back(this->GetValue(i) - control_point.GetValue(i));
  }
  return ControlPoint(coordinates_new);
}

baf::ControlPoint baf::ControlPoint::operator*(const double &scalar) const {
  std::vector<double> coordinates_new(this->GetDimension());
  std::transform(coordinates_.begin(), coordinates_.end(), coordinates_new.begin(),
      std::bind(std::multiplies<>(), std::placeholders::_1, scalar));
  return ControlPoint(coordinates_new);
}

baf::ControlPoint baf::ControlPoint::Transform(std::array<std::array<double, 4>, 4> TransMatrix,
    std::array<double, 3> scaling) const {
  std::vector<double> coordinates_new = {TransMatrix[0][3], TransMatrix[1][3], TransMatrix[2][3]};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      coordinates_new[j] += TransMatrix[j][i] * scaling[i] * GetValue(i);
    }
  }
  return ControlPoint(coordinates_new);
}

double baf::ControlPoint::GetEuclideanNorm() const {
  double euclidean_norm = 0.0;
  for (auto &coordinate : coordinates_) {
    euclidean_norm += pow(coordinate, 2);
  }
  euclidean_norm = sqrt(euclidean_norm);
  return euclidean_norm;
}
