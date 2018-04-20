/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "integration_rule_5_points.h"

IntegrationRule5Points::IntegrationRule5Points() :
    OneDimensionalIntegrationRule({-(1.0 / 3) * sqrt(5 + 2.0 * sqrt(10.0 / 7)),
                                   -(1.0 / 3) * sqrt(5 - 2.0 * sqrt(10.0 / 7)), 0,
                                   (1.0 / 3) * sqrt(5 - 2.0 * sqrt(10.0 / 7)),
                                   (1.0 / 3) * sqrt(5 + 2.0 * sqrt(10.0 / 7))},
                                  {(322.0 - 13 * sqrt(70)) / 900, (322.0 + 13 * sqrt(70)) / 900, 128.0 / 225,
                                   (322.0 + 13 * sqrt(70)) / 900, (322.0 - 13 * sqrt(70)) / 900}) {}
