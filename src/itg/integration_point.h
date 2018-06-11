/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <array>

#ifndef SRC_ITG_INTEGRATION_POINT_H_
#define SRC_ITG_INTEGRATION_POINT_H_

namespace itg {
template<int DIM>
class IntegrationPoint {
 public:
  IntegrationPoint(const std::array<double, DIM> &coordinates, double weight)
      : coordinates_(coordinates), weight_(weight) {}

  int GetDimension() const {
    return DIM;
  }

  std::array<double, DIM> GetCoordinates() const {
    return coordinates_;
  }

  double GetWeight() const {
    return weight_;
  }

 private:
  std::array<double, DIM> coordinates_;
  double weight_;
};
}  // namespace itg

#endif  // SRC_ITG_INTEGRATION_POINT_H_
