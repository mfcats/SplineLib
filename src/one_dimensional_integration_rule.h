/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_ONE_DIMENSIONAL_INTEGRATION_RULE_H
#define SPLINELIB_ONE_DIMENSIONAL_INTEGRATION_RULE_H

#include <stdexcept>
#include <vector>
#include <cmath>

class OneDimensionalIntegrationRule {
 public:
  explicit OneDimensionalIntegrationRule(int points);

  int points();
  double point(int point);
  double weight(int point);

 private:
  int number_of_points_;
  std::vector<double> points_;
  std::vector<double> weights_;
};

#endif //SPLINELIB_ONE_DIMENSIONAL_INTEGRATION_RULE_H
