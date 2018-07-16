/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_ITG_FOUR_POINT_GAUSS_LEGENDRE_H_
#define SRC_ITG_FOUR_POINT_GAUSS_LEGENDRE_H_

#include <cmath>
#include <array>

#include "integration_rule.h"

namespace itg {
template<int DIM>
class FourPointGaussLegendre : public IntegrationRule<DIM> {
 public:
  FourPointGaussLegendre() : IntegrationRule<DIM>({IntegrationPoint<1>({points_[0]}, weights_[0]),
                                                   IntegrationPoint<1>({points_[1]}, weights_[1]),
                                                   IntegrationPoint<1>({points_[2]}, weights_[2]),
                                                   IntegrationPoint<1>({points_[3]}, weights_[3])}) {}

 private:
  static constexpr std::array<double, 4> weights_ =
      {(18.0 - sqrt(30)) / 36, (18.0 + sqrt(30)) / 36, (18.0 + sqrt(30)) / 36, (18.0 - sqrt(30)) / 36};
  static constexpr std::array<double, 4> points_ =
      {-sqrt(3.0 / 7 + 2.0 / 7 * sqrt(6.0 / 5)), -sqrt(3.0 / 7 - 2.0 / 7 * sqrt(6.0 / 5)),
       sqrt(3.0 / 7 - 2.0 / 7 * sqrt(6.0 / 5)), sqrt(3.0 / 7 + 2.0 / 7 * sqrt(6.0 / 5))};
};
}  // namespace itg

#endif  // SRC_ITG_FOUR_POINT_GAUSS_LEGENDRE_H_
