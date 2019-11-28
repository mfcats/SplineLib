/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

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

class AnIntegrationRule : public Test {
 public:
  AnIntegrationRule() {
    rules_.emplace_back(iga::itg::OnePointGaussLegendre());
    rules_.emplace_back(iga::itg::TwoPointGaussLegendre());
    rules_.emplace_back(iga::itg::ThreePointGaussLegendre());
    rules_.emplace_back(iga::itg::FourPointGaussLegendre());
    rules_.emplace_back(iga::itg::FivePointGaussLegendre());
  }

 protected:
  std::vector<iga::itg::IntegrationRule> rules_;
};

TEST_F(AnIntegrationRule, ReturnsCorrectNumberOfPoints) {  // NOLINT
  for (int points = 1; points <= 5; points++) {
    ASSERT_THAT(rules_[points - 1].GetNumberOfIntegrationPoints(), points);
  }
}

TEST_F(AnIntegrationRule, ReturnsCorrectWeightSum) {  // NOLINT
  for (int points = 1; points <= 5; points++) {
    double weight_sum = 0;
    for (int weight = 0; weight < points; weight++) {
      weight_sum += rules_[points - 1].GetIntegrationPoints()[weight].GetWeight();
    }
    ASSERT_THAT(weight_sum, DoubleEq(2));
  }
}

TEST_F(AnIntegrationRule, ReturnsCorrectPointSum) {  // NOLINT
  for (int points = 1; points <= 5; points++) {
    double point_sum = 0;
    for (int point = 0; point < points; point++) {
      point_sum += rules_[points - 1].GetIntegrationPoints()[point].GetCoordinate();
    }
    ASSERT_THAT(point_sum, DoubleNear(0.0, util::numeric_settings::GetEpsilon<double>()));
  }
}
