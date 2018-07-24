/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_ITG_INTEGRATION_RULE_H_
#define SRC_ITG_INTEGRATION_RULE_H_

#include <cmath>
#include <vector>

#include "integration_point.h"
#include "multi_index_handler.h"

namespace itg {
template<int DIM>
class IntegrationRule {
 public:
  virtual ~IntegrationRule() = default;

  explicit IntegrationRule(const std::vector<IntegrationPoint<1>> &points) : points_(points) {}

  int GetNumberOfIntegrationPoints() const {
    return pow(points_.size(), DIM);
  }

  double GetCoordinate(int point, int dimension) const {
#ifdef DEBUG
    return points_.at(point).GetCoordinates().at(dimension);
#else
    return points_[point].GetCoordinates()[dimension];
#endif
  }

  std::vector<IntegrationPoint<DIM>> GetIntegrationPoints() {
    std::vector<IntegrationPoint<DIM>> integration_points;
    std::array<int, DIM> max_dimension_points;
    for (int i = 0; i < DIM; i++) {
      max_dimension_points[i] = points_.size();
    }
    util::MultiIndexHandler<DIM> multiIndexHandler(max_dimension_points);
    for (int i = 0; i < GetNumberOfIntegrationPoints(); i++) {
      double weight = 1;
      std::array<double, DIM> coordinates = {0};
      for (int j = 0; j < DIM; j++) {
        weight *= points_[multiIndexHandler[j]].GetWeight();
        coordinates[j] = points_[multiIndexHandler[j]].GetCoordinates()[j];
      }
      ++multiIndexHandler;
      integration_points.push_back(IntegrationPoint<DIM>(coordinates, weight));
    }
    return integration_points;
  }

 private:
  std::vector<IntegrationPoint<1>> points_;
};
}  // namespace itg

#endif  // SRC_ITG_INTEGRATION_RULE_H_
