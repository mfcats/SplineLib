/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_BASIS_FUNCTION_H_
#define SRC_IGA_BASIS_FUNCTION_H_

#include "element_integration_point_iga.h"
#include "spline.h"

namespace iga {
class BasisFunctionHandler {
 public:
  explicit BasisFunctionHandler(std::shared_ptr<spl::Spline<2>> spl) : spline_(std::move(spl)) {}

  std::vector<iga::ElementIntegrationPoint>
  EvaluateAllElementNonZeroBasisFunctions(int direction,
                                          int element_number,
                                          const iga::IntegrationRule<1> &rule) const {
    iga::Element element = GetElementList(direction)[element_number];
    std::vector<iga::ElementIntegrationPoint> element_integration_points;
    std::vector<double> non_zero_basis_functions;
    for (int i = 0; i < rule.GetNumberOfIntegrationPoints(); ++i) {
      ParamCoord integration_point =
          ReferenceSpace2ParameterSpace(element.GetNode(1), element.GetNode(0), rule.GetCoordinate(i, 0));
      non_zero_basis_functions = EvaluateAllNonZeroBasisFunctions(direction, ParamCoord{integration_point});
      element_integration_points.emplace_back(iga::ElementIntegrationPoint(non_zero_basis_functions));
    }
    return element_integration_points;
  }

 private:
  std::shared_ptr<spl::Spline<2>> spline_;
};
}

#endif  // SRC_IGA_BASIS_FUNCTION_H_
