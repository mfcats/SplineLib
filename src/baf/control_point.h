/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_BAF_CONTROL_POINT_H_
#define SRC_BAF_CONTROL_POINT_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <vector>

namespace baf {
class ControlPoint {
 public:
  explicit ControlPoint(std::initializer_list<double> coordinates);
  explicit ControlPoint(std::vector<double> coordinates);
  explicit ControlPoint(uint64_t dimension);

  int GetDimension() const;
  double GetValue(int dimension) const;

  void SetValue(int dimension, double value);

  ControlPoint operator+(const ControlPoint &control_point) const;
  ControlPoint operator-(const ControlPoint &control_point) const;
  ControlPoint operator*(const double &scalar) const;

  ControlPoint Transform(std::array<std::array<double, 4>, 4> TransMatrix,
                                       std::array<double, 3> scaling) const;
  double GetEuclideanNorm() const;

 protected:
  std::vector<double> coordinates_;
};
}  // namespace baf

#endif  // SRC_BAF_CONTROL_POINT_H_
