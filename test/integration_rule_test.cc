/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "integration_rule.h"
#include <gmock/gmock.h>

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class A1DIntegrationRule : public Test {
 public:
  A1DIntegrationRule() {
    for (int points = 1; points <= 5; points++) {
      rules_.emplace_back(OneDimensionalIntegrationRule(points));
    }
  };

 protected:
  std::vector<OneDimensionalIntegrationRule> rules_;
};

TEST_F(A1DIntegrationRule, ReturnsCorrectNumberOfPoints) {
  for (int points = 1; points <= 5; points++) {
    ASSERT_THAT(rules_[points - 1].points(), points);
  }
}

TEST_F(A1DIntegrationRule, ReturnsCorrectWeights) {
  for (int points = 1; points <= 5; points++) {
    double weight_sum = 0;
    for (int weight = 0; weight < points; weight++) {
      weight_sum += rules_[points - 1].weight(weight);
    }
    ASSERT_THAT(weight_sum, DoubleEq(2));
  }
}

TEST_F(A1DIntegrationRule, ReturnsCorrectPoints) {
  for (int points = 1; points <= 3; points++) {
    double point_sum = 0;
    for (int point = 0; point < points; point++) {
      point_sum += rules_[points - 1].point(point);
    }
    ASSERT_THAT(point_sum, DoubleNear(0, pow(10, -35)));
  }
}

class A2DIntegrationRuleWith3PointsEach : public Test {
 public:
  A2DIntegrationRuleWith3PointsEach() : rule_(IntegrationRule<2>(3)) {}

 protected:
  IntegrationRule<2> rule_;
};

TEST_F(A2DIntegrationRuleWith3PointsEach, ReturnsCorrectNumberOfPoints) {
  ASSERT_THAT(rule_.points(), 9);
}

TEST_F(A2DIntegrationRuleWith3PointsEach, ReturnsCorrectPoint) {
  ASSERT_THAT(rule_.point(1, 0), 0);
  ASSERT_THAT(rule_.point(1, 1), 0);
}

TEST_F(A2DIntegrationRuleWith3PointsEach, ReturnsCorrectWeights) {
  ASSERT_THAT(rule_.weight(1, 0), DoubleEq(8.0 / 9.0));
  ASSERT_THAT(rule_.weight(0, 1), DoubleEq(5.0 / 9.0));
}

TEST(A3DIntegrationRuleWith3PointsEach, ReturnsCorrectNumberOfPoints) {
  ASSERT_THAT(IntegrationRule<3>(3).points(), 27);
}

TEST(A3DIntegrationRuleWith1PointEach, ReturnsCorrectNumberOfPoints) {
  ASSERT_THAT(IntegrationRule<3>(1).points(), 1);
}
