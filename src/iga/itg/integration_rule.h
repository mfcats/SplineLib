/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_ITG_INTEGRATION_RULE_H_
#define SRC_IGA_ITG_INTEGRATION_RULE_H_

#include <cmath>
#include <vector>

#include "integration_point.h"
#include "multi_index_handler.h"

namespace iga {
namespace itg {
class IntegrationRule {
 public:
  virtual ~IntegrationRule() = default;

  explicit IntegrationRule(std::vector<IntegrationPoint> points) : points_(std::move(points)) {}

  int GetNumberOfIntegrationPoints() const {
    return points_.size();
  }

  double GetCoordinate(int point) const {
#ifdef DEBUG
    return points_.at(point).GetCoordinate();
#else
    return points_.at(point).GetCoordinate();
#endif
  }

  std::vector<IntegrationPoint> GetIntegrationPoints() {
    return points_;
  }

 private:
  std::vector<IntegrationPoint> points_;
};
}  // namespace itg
}  // namespace iga

#endif  // SRC_IGA_ITG_INTEGRATION_RULE_H_
