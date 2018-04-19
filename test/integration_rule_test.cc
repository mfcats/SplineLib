/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "one_dimensional_integration_rule.h"

#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::DoubleNear;

TEST(A1DIntegrationRule, ReturnsCorrectNumberOfPoints) {
  for (int points = 1; points <= 5; points++) {
    ASSERT_THAT(OneDimensionalIntegrationRule(points).points(), points);
  }
}

TEST(A1DIntegrationRule, ReturnsCorrectWeights) {
  for (int points = 1; points <= 5; points++) {
    OneDimensionalIntegrationRule rule = OneDimensionalIntegrationRule(points);
    double weight_sum = 0;
    for (int weight = 0; weight < points; weight++) {
      weight_sum += rule.weight(weight);
    }
    ASSERT_THAT(weight_sum, DoubleEq(2));
  }
}

TEST(A1DIntegrationRule, ReturnsCorrectPoints) {
  for (int points = 1; points <= 3; points++) {
    OneDimensionalIntegrationRule rule = OneDimensionalIntegrationRule(points);
    double point_sum = 0;
    for (int point = 0; point < points; point++) {
      point_sum += rule.point(point);
    }
    ASSERT_THAT(point_sum, DoubleNear(0, pow(10, -35)));
  }
}
