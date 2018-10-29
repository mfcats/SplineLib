/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_FOUR_POINT_GAUSS_LEGENDRE_H_
#define SRC_IGA_FOUR_POINT_GAUSS_LEGENDRE_H_

#include <cmath>
#include <array>

#include "integration_rule.h"

namespace iga {
namespace itg {
template<int DIM>
class FourPointGaussLegendre : public IntegrationRule<DIM> {
 public:
  FourPointGaussLegendre() : IntegrationRule<DIM>({IntegrationPoint<1>(std::array<double, 1>{
      -sqrt(3.0 / 7 + 2.0 / 7 * sqrt(6.0 / 5))}, (18.0 - sqrt(30)) / 36),
                                                   IntegrationPoint<1>(std::array<double, 1>{
                                                                           -sqrt(3.0 / 7 - 2.0 / 7 * sqrt(6.0 / 5))},
                                                                       (18.0 + sqrt(30)) / 36),
                                                   IntegrationPoint<1>(std::array<double, 1>{
                                                                           sqrt(3.0 / 7 - 2.0 / 7 * sqrt(6.0 / 5))},
                                                                       (18.0 + sqrt(30)) / 36),
                                                   IntegrationPoint<1>(std::array<double, 1>{
                                                                           sqrt(3.0 / 7 + 2.0 / 7 * sqrt(6.0 / 5))},
                                                                       (18.0 - sqrt(30)) / 36)}) {}
};
}  // namespace itg
}  // namespace iga

#endif  // SRC_IGA_FOUR_POINT_GAUSS_LEGENDRE_H_
