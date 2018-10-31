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

#ifndef SRC_IGA_ITG_INTEGRATION_POINT_H_
#define SRC_IGA_ITG_INTEGRATION_POINT_H_

namespace iga {
namespace itg {
class IntegrationPoint {
 public:
  IntegrationPoint(double coordinate, double weight)
      : coordinate_(coordinate), weight_(weight) {}

  double GetCoordinate() const {
    return coordinate_;
  }

  double GetWeight() const {
    return weight_;
  }

 private:
  double coordinate_;
  double weight_;
};
}  // namespace itg
}  // namespace iga

#endif  // SRC_IGA_ITG_INTEGRATION_POINT_H_
