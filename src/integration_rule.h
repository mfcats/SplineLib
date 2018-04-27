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

#include "one_dimensional_integration_rule.h"

template<int dimensions>
class IntegrationRule {
 public:
  explicit IntegrationRule(const OneDimensionalIntegrationRule &rule) : rules_(rule) {}

  int points() const {
    return pow(rules_.points(), dimensions);
  }

  double point(int point, int dimension) const {
    return rules_.point(point);
  }

  double weight(int point, int dimension) const {
    return rules_.weight(point);
  }

 private:
  OneDimensionalIntegrationRule rules_;
};

#endif  // SRC_INTEGRATION_RULE_H_
