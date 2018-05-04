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

#include "five_point_gauss_legendre.h"
#include "four_point_gauss_legendre.h"
#include "integration_rule.h"
#include "numeric_settings.h"
#include "one_point_gauss_legendre.h"
#include "three_point_gauss_legendre.h"
#include "two_point_gauss_legendre.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class A1DIntegrationRule : public Test {
 public:
  A1DIntegrationRule() {
    rules_.emplace_back(OnePointGaussLegendre<1>());
    rules_.emplace_back(TwoPointGaussLegendre<1>());
    rules_.emplace_back(ThreePointGaussLegendre<1>());
    rules_.emplace_back(FourPointGaussLegendre<1>());
    rules_.emplace_back(FivePointGaussLegendre<1>());
  };

 protected:
  std::vector<IntegrationRule<1>> rules_;
};

TEST_F(A1DIntegrationRule, ReturnsCorrectNumberOfPoints) {
  for (int points = 1; points <= 5; points++) {
    ASSERT_THAT(rules_[points - 1].GetNumberOfIntegrationPoints(), points);
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
  A2DIntegrationRuleWith3Points() : rule_(ThreePointGaussLegendre<2>()) {}

 protected:
  IntegrationRule<2> rule_;
};

TEST_F(A2DIntegrationRuleWith3Points, ReturnsCorrectNumberOfPoints) {
  ASSERT_THAT(rule_.GetNumberOfIntegrationPoints(), 9);
  ASSERT_THAT(rule_.GetIntegrationPoints().size(), 9);
}

TEST_F(A2DIntegrationRuleWith3Points, ReturnsCorrectPoint) {
  ASSERT_THAT(rule_.coordinate(0, 0), DoubleNear(-sqrt(3.0 / 5), NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(rule_.coordinate(1, 0), DoubleNear(0, NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(rule_.coordinate(2, 0), DoubleNear(sqrt(3.0 / 5), NumericSettings<double>::kEpsilon()));
}

TEST_F(A2DIntegrationRuleWith3Points, ReturnsCorrectWeightSum) {
  double weight_sum = 0;
  for (int weight = 0; weight < rule_.GetNumberOfIntegrationPoints(); weight++) {
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
  ASSERT_THAT(ThreePointGaussLegendre<3>().GetNumberOfIntegrationPoints(), 27);
}

TEST(A3DIntegrationRuleWith1Point, ReturnsCorrectNumberOfPoints) {
  ASSERT_THAT(OnePointGaussLegendre<3>().GetNumberOfIntegrationPoints(), 1);
}
