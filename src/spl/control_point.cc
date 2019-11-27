/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "src/spl/control_point.h"

#include <cmath>
#include <functional>
#include <utility>

#include "src/util/stl_container_access.h"
#include "src/util/vector_utils.h"

namespace splinelib::src::baf {
ControlPoint::ControlPoint(std::vector<double> coordinates) : coordinates_(std::move(coordinates)) {}

ControlPoint::ControlPoint(std::initializer_list<double> const &coordinates) : coordinates_(coordinates) {}

ControlPoint::ControlPoint(int number_of_dimensions) : coordinates_(std::vector(number_of_dimensions, 0.0)) {}

ControlPoint::ControlPoint(ControlPoint &&other) noexcept : coordinates_(std::move(other.coordinates_)) {}

ControlPoint & ControlPoint::operator=(ControlPoint &&rhs) noexcept {
  coordinates_ = std::move(rhs.coordinates_);
  return (*this);
}

ControlPoint ControlPoint::Transform(std::array<std::array<double, 4>, 4> const &transformation_matrix,
                                     std::array<double, 3> const &scaling) const {
  if (coordinates_.size() < 3) throw std::logic_error("splinelib::src::spl::ControlPoint::Transform: only vectors with"
                                                       "dimension of at least three can be transformed.");
  std::vector<double> coordinates_new = {GetValue(GetValue(transformation_matrix, 0), 3),
                                         GetValue(GetValue(transformation_matrix, 1), 3),
                                         GetValue(GetValue(transformation_matrix, 2), 3)};
  for (Dimension i{0}; i < Dimension(3); ++i) {
    for (Dimension j{0}; j < Dimension{3}; ++j) {
      GetValue(coordinates_new, j) += GetValue(GetValue(transformation_matrix, j), i) * GetValue(scaling, i) *
                                      GetValue(coordinates_, i);
    }
  }
  return ControlPoint(coordinates_new);
}

double ControlPoint::GetEuclideanNorm() const {
  double euclidean_norm = 0.0;
  for (auto &coordinate : coordinates_) {
    euclidean_norm += pow(coordinate, 2);
  }
  euclidean_norm = sqrt(euclidean_norm);
  return euclidean_norm;
}

ControlPoint operator+(ControlPoint const &lhs, ControlPoint const &rhs) {
  std::vector<double> const coordinates_new = util::vector_utils::ComputeSum(lhs.coordinates_, rhs.coordinates_);
  return ControlPoint(coordinates_new);
}

ControlPoint operator-(ControlPoint const &lhs, ControlPoint const &rhs) {
  std::vector<double> const coordinates_new = util::vector_utils::ComputeDifference(lhs.coordinates_, rhs.coordinates_);
  return ControlPoint(coordinates_new);
}

ControlPoint operator*(ControlPoint const &control_point, double scalar) {
  std::vector<double> coordinates_new(control_point.GetDimensionality());
  std::transform(control_point.coordinates_.begin(), control_point.coordinates_.end(), coordinates_new.begin(),
                 std::bind(std::multiplies<>(), std::placeholders::_1, scalar));
  return ControlPoint(coordinates_new);
}

ControlPoint operator*(double scalar, ControlPoint const &control_point) {
  std::vector<double> coordinates_new(control_point.GetDimensionality());
  std::transform(control_point.coordinates_.begin(), control_point.coordinates_.end(), coordinates_new.begin(),
                 std::bind(std::multiplies<>(), std::placeholders::_1, scalar));
  return ControlPoint(coordinates_new);
}
}  // namespace splinelib::src::baf
