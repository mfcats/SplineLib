/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_INTEGRATION_RULE_H_
#define SRC_INTEGRATION_RULE_H_

#include <vector>
#include <cmath>

#include "integration_point.h"
#include "multi_index_handler.h"

template<int dimensions>
class IntegrationRule {
 public:
  explicit IntegrationRule(const std::vector<IntegrationPoint<1>> &points) : points_(points) {}

  int points() const {
    return pow(points_.size(), dimensions);
  }

  int point(int point, int dimension) const {
    return points_[0].GetCoordinates()[dimension];
  }

  std::vector<IntegrationPoint<dimensions>> GetIntegrationPoints() {
    std::vector<IntegrationPoint<dimensions>> integration_points;
    std::array<int, dimensions> max_dimension_points;
    for (int i = 0; i < dimensions; i++) {
      max_dimension_points[i] = dimensions;
    }
    MultiIndexHandler<dimensions> multiIndexHandler(max_dimension_points);
    for (int i = 0; i < points(); i++) {
      double weight = 1;
      std::array<double, dimensions> coordinates;
      for (int j = 0; j < dimensions; j++) {
        weight *= points_[multiIndexHandler[j]].GetWeight();
        coordinates[j] = points_[multiIndexHandler[j]].GetCoordinates()[j];
      }
      multiIndexHandler++;
      integration_points.push_back(IntegrationPoint<dimensions>(coordinates, weight));
    }
    return integration_points;
  }

 private:
  std::vector<IntegrationPoint<1>> points_;
};

#endif  // SRC_INTEGRATION_RULE_H_
