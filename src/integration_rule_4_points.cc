/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "integration_rule_4_points.h"

IntegrationRule4Points::IntegrationRule4Points() :
    OneDimensionalIntegrationRule({-sqrt(3.0 / 7 + 2.0 / 7 * sqrt(6.0 / 5)), -sqrt(3.0 / 7 - 2.0 / 7 * sqrt(6.0 / 5)),
                                   sqrt(3.0 / 7 - 2.0 / 7 * sqrt(6.0 / 5)), sqrt(3.0 / 7 + 2.0 / 7 * sqrt(6.0 / 5))},
                                  {(18.0 - sqrt(30)) / 36, (18.0 + sqrt(30)) / 36, (18.0 + sqrt(30)) / 36,
                                   (18.0 - sqrt(30)) / 36}) {}
