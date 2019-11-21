/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_SPL_CONTROL_POINT_H_
#define SRC_SPL_CONTROL_POINT_H_

#include <algorithm>
#include <array>
#include <functional>
#include <initializer_list>
#include <vector>

#include "src/util/named_type.h"

namespace splinelib::src::baf {
class ControlPoint {
 public:
  explicit ControlPoint(std::vector<double> coordinates);
  ControlPoint(std::initializer_list<double> const &coordinates);
  explicit ControlPoint(uint64_t dimension);

  ControlPoint(ControlPoint const &other) = default;
  ControlPoint(ControlPoint &&other) noexcept;
  ControlPoint & operator=(ControlPoint const &rhs) = default;
  ControlPoint & operator=(ControlPoint &&rhs) noexcept;

  int GetDimensionality() const;
  double GetValueForDimension(Dimension const &dimension) const;
  double operator[](Dimension const &dimension) const;

  void SetValue(int dimension, double value);

  friend ControlPoint operator+(ControlPoint const &lhs, ControlPoint const &rhs);
  friend ControlPoint operator-(ControlPoint const &lhs, ControlPoint const &rhs);
  friend ControlPoint operator*(ControlPoint const &control_point, double scalar);
  friend ControlPoint operator*(double scalar, ControlPoint const &control_point);

  ControlPoint Transform(std::array<std::array<double, 4>, 4> const &TransformationMatrix,
                         std::array<double, 3> const &scaling) const;

  double GetEuclideanNorm() const;

 protected:
  std::vector<double> coordinates_;
};

ControlPoint operator+(ControlPoint const &lhs, ControlPoint const &rhs);
ControlPoint operator-(ControlPoint const &lhs, ControlPoint const &rhs);
ControlPoint operator*(ControlPoint const &control_point, double scalar);
ControlPoint operator*(double scalar, ControlPoint const &control_point);

#include "src/spl/control_point.inc"
}  // namespace splinelib::src::baf

#endif  // SRC_SPL_CONTROL_POINT_H_
