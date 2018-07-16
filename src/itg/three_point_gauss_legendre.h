/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_ITG_THREE_POINT_GAUSS_LEGENDRE_H_
#define SRC_ITG_THREE_POINT_GAUSS_LEGENDRE_H_

#include <cmath>
#include <array>

#include "integration_rule.h"

namespace itg {
template<int DIM>
class ThreePointGaussLegendre : public IntegrationRule<DIM> {
 public:
  ThreePointGaussLegendre() : IntegrationRule<DIM>({IntegrationPoint<1>({points_[0]}, weights_[0]),
                                                    IntegrationPoint<1>({points_[1]}, weights_[1]),
                                                    IntegrationPoint<1>({points_[2]}, weights_[2])}) {}

 private:
  static constexpr std::array<double, 3> weights_ = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
  static constexpr std::array<double, 3> points_ = {-sqrt(3.0 / 5), 0, sqrt(3.0 / 5)};
};
}  // namespace itg

#endif  // SRC_ITG_THREE_POINT_GAUSS_LEGENDRE_H_
