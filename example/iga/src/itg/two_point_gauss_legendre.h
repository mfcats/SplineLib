/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_IGA_ITG_TWO_POINT_GAUSS_LEGENDRE_H_
#define SRC_IGA_ITG_TWO_POINT_GAUSS_LEGENDRE_H_

#include <cmath>
#include <array>

#include "integration_rule.h"

namespace iga {
namespace itg {
class TwoPointGaussLegendre : public IntegrationRule {
 public:
  TwoPointGaussLegendre() : IntegrationRule({IntegrationPoint(-sqrt(1.0 / 3), 1),
                                             IntegrationPoint(sqrt(1.0 / 3), 1)}) {}
};
}  // namespace itg
}  // namespace iga

#endif  // SRC_IGA_ITG_TWO_POINT_GAUSS_LEGENDRE_H_
