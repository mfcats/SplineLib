/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "gmock/gmock.h"

#include <cmath>

#include "integration_rule.h"
#include "integration_rule_1_point.h"
#include "integration_rule_2_points.h"
#include "integration_rule_3_points.h"
#include "integration_rule_4_points.h"
#include "integration_rule_5_points.h"
#include "numeric_settings.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class A1DIntegrationPoint : public Test {
 public:
  A1DIntegrationPoint() : integration_point_({1.5}, 0.5) {}

 protected:
  IntegrationPoint<1> integration_point_;
};

TEST_F(A1DIntegrationPoint, ReturnsCorrectCoordinate) {
  ASSERT_THAT(integration_point_.GetCoordinates().size(), 1);
  ASSERT_THAT(integration_point_.GetCoordinates()[0], DoubleEq(1.5));
}

TEST_F(A1DIntegrationPoint, ReturnsCorrectWeight) {
  ASSERT_THAT(integration_point_.GetWeight(), DoubleEq(0.5));
}

TEST_F(A1DIntegrationPoint, ReturnsCorrectDimension) {
  ASSERT_THAT(integration_point_.GetDimension(), 1);
}

class A1DIntegrationRule : public Test {
 public:
  A1DIntegrationRule() {
    rules_.emplace_back(IntegrationRule1Point<1>());
    rules_.emplace_back(IntegrationRule2Points<1>());
    rules_.emplace_back(IntegrationRule3Points<1>());
    rules_.emplace_back(IntegrationRule4Points<1>());
    rules_.emplace_back(IntegrationRule5Points<1>());
  };

 protected:
  std::vector<IntegrationRule<1>> rules_;
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
      weight_sum += rules_[points - 1].GetIntegrationPoints()[weight].GetWeight();
    }
    //ASSERT_THAT((2 * (322.0 - 13 * sqrt(70)) / 900.0 + 2*(322.0 + 13 * sqrt(70)) / 900.0 + 128.0 / 225.0), DoubleEq(2));
    ASSERT_THAT(weight_sum, DoubleEq(2));
  }
}

TEST_F(A1DIntegrationRule, ReturnsCorrectPoints) {
  for (int points = 1; points <= 5; points++) {
    double point_sum = 0;
    for (int point = 0; point < points; point++) {
      point_sum += rules_[points - 1].GetIntegrationPoints()[point].GetCoordinates()[0];
    }
    ASSERT_THAT(point_sum, DoubleNear(0.0, NumericSettings<double>::kEpsilon()));
  }
}

class A2DIntegrationRuleWith3Points : public Test {
 public:
  A2DIntegrationRuleWith3Points() : rule_(IntegrationRule3Points<2>()) {}

 protected:
  IntegrationRule<2> rule_;
};

TEST_F(A2DIntegrationRuleWith3Points, ReturnsCorrectNumberOfPoints) {
  ASSERT_THAT(rule_.points(), 9);
}

TEST_F(A2DIntegrationRuleWith3Points, ReturnsCorrectPoint) {
  ASSERT_THAT(rule_.point(1, 0), 0);
  ASSERT_THAT(rule_.point(1, 1), 0);
}

TEST_F(A2DIntegrationRuleWith3Points, ReturnsCorrectWeights) {
  ASSERT_THAT(rule_.GetIntegrationPoints()[0].GetWeight(), DoubleEq(25.0 / 81.0));
  ASSERT_THAT(rule_.GetIntegrationPoints()[1].GetWeight(), DoubleEq(40.0 / 81.0));
  ASSERT_THAT(rule_.GetIntegrationPoints()[2].GetWeight(), DoubleEq(25.0 / 81.0));
  ASSERT_THAT(rule_.GetIntegrationPoints()[3].GetWeight(), DoubleEq(40.0 / 81.0));
  ASSERT_THAT(rule_.GetIntegrationPoints()[4].GetWeight(), DoubleEq(64.0 / 81.0));
  ASSERT_THAT(rule_.GetIntegrationPoints()[5].GetWeight(), DoubleEq(40.0 / 81.0));
  ASSERT_THAT(rule_.GetIntegrationPoints()[6].GetWeight(), DoubleEq(25.0 / 81.0));
  ASSERT_THAT(rule_.GetIntegrationPoints()[7].GetWeight(), DoubleEq(40.0 / 81.0));
  ASSERT_THAT(rule_.GetIntegrationPoints()[8].GetWeight(), DoubleEq(25.0 / 81.0));
}

TEST(A3DIntegrationRuleWith3Points, ReturnsCorrectNumberOfPoints) {
  ASSERT_THAT(IntegrationRule3Points<3>().points(), 27);
}

TEST(A3DIntegrationRuleWith1Point, ReturnsCorrectNumberOfPoints) {
  ASSERT_THAT(IntegrationRule1Point<3>().points(), 1);
}
