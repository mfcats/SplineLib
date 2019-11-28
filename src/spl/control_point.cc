/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "src/spl/control_point.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>

#include "src/util/vector_utils.h"

namespace splinelib::src::spl {
ControlPoint::ControlPoint(std::vector<double> coordinates) : coordinates_(std::move(coordinates)) {}

ControlPoint::ControlPoint(std::initializer_list<double> const &coordinates) : coordinates_(coordinates) {}

ControlPoint::ControlPoint(int number_of_dimensions) : coordinates_(std::vector(number_of_dimensions, 0.0)) {}

ControlPoint::ControlPoint(ControlPoint &&other) noexcept : coordinates_(std::move(other.coordinates_)) {}

ControlPoint & ControlPoint::operator=(ControlPoint &&rhs) noexcept {
  coordinates_ = std::move(rhs.coordinates_);
  return (*this);
}

double ControlPoint::ComputeTwoNorm() const {
  return util::vector_utils::ComputeTwoNorm(coordinates_);
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
}  // namespace splinelib::src::spl
