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

TEST_F(A1DIntegrationRule, ReturnsCorrectWeightSum) {
  for (int points = 1; points <= 5; points++) {
    double weight_sum = 0;
    for (int weight = 0; weight < points; weight++) {
      weight_sum += rules_[points - 1].GetIntegrationPoints()[weight].GetWeight();
    }
    ASSERT_THAT(weight_sum, DoubleEq(2));
  }
}

TEST_F(A1DIntegrationRule, ReturnsCorrectPointSum) {
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
  ASSERT_THAT(rule_.GetIntegrationPoints().size(), 9);
}

TEST_F(A2DIntegrationRuleWith3Points, ReturnsCorrectPoint) {
  ASSERT_THAT(rule_.coordinate(0, 0), DoubleNear(-sqrt(3.0 / 5), NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(rule_.coordinate(1, 0), DoubleNear(0, NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(rule_.coordinate(2, 0), DoubleNear(sqrt(3.0 / 5), NumericSettings<double>::kEpsilon()));
}

TEST_F(A2DIntegrationRuleWith3Points, ReturnsCorrectWeightSum) {
  double weight_sum = 0;
  for (int weight = 0; weight < rule_.points(); weight++) {
    weight_sum += rule_.GetIntegrationPoints()[weight].GetWeight();
  }
  ASSERT_THAT(weight_sum, DoubleEq(4));
}

TEST_F(A2DIntegrationRuleWith3Points, ReturnsCorrectPointSum) {
  double point_sum = 0;
  for (int weight_dim0 = 0; weight_dim0 < 3; weight_dim0++) {
    for (int weight_dim1 = 0; weight_dim1 < 3; weight_dim1++) {
      for (int dimension = 0; dimension < 2; dimension++) {
        point_sum += rule_.GetIntegrationPoints()[weight_dim1 * 3 + weight_dim0].GetCoordinates()[0];
      }
    }
  }
  ASSERT_THAT(point_sum, DoubleEq(0));
}

TEST(A3DIntegrationRuleWith3Points, ReturnsCorrectNumberOfPoints) {
  ASSERT_THAT(IntegrationRule3Points<3>().points(), 27);
}

TEST(A3DIntegrationRuleWith1Point, ReturnsCorrectNumberOfPoints) {
  ASSERT_THAT(IntegrationRule1Point<3>().points(), 1);
}
