/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_BASIS_FUNCTION_HANDLER_H_
#define SRC_IGA_BASIS_FUNCTION_HANDLER_H_

#include <math.h>
#include <array>
#include <vector>

#include "connectivity_handler.h"
#include "element.h"
#include "element_generator.h"
#include "element_integration_point.h"
#include "integration_rule.h"
#include "mapping_handler.h"
#include "nurbs.h"

namespace iga {
class BasisFunctionHandler {
 public:
  explicit BasisFunctionHandler(std::shared_ptr<spl::NURBS<2>> spl);

  std::vector<iga::elm::ElementIntegrationPoint>
  EvaluateAllElementNonZeroNURBSBasisFunctions(int element_number, const iga::itg::IntegrationRule &rule) const;

  std::vector<iga::elm::ElementIntegrationPoint> EvaluateAllElementNonZeroNURBSBasisFunctionDerivatives(
      int element_number, const iga::itg::IntegrationRule &rule) const;

  std::vector<iga::elm::ElementIntegrationPoint> EvaluateAllElementNonZeroNURBSBafDerivativesPhysical(
      int element_number, const iga::itg::IntegrationRule &rule) const;

 private:
  std::vector<double> EvaluateAllNonZeroNURBSBasisFunctions(std::array<ParamCoord, 2> param_coord) const;

  std::array<std::vector<double>, 2> EvaluateAllNonZeroNURBSBasisFunctionDerivatives(
      std::array<ParamCoord, 2> param_coord) const;

  std::array<std::vector<double>, 2> EvaluateAllNonZeroNURBSBafDerivativesPhyiscal(
      std::array<ParamCoord, 2> param_coord) const;

  double GetWeight(std::array<ParamCoord, 2> param_coord, int local_index) const;

  std::shared_ptr<spl::NURBS<2>> spline_;
  std::shared_ptr<iga::MappingHandler> mapping_handler_;
  std::shared_ptr<iga::elm::ElementGenerator<2>> element_generator_;
};
}  // namespace iga

#endif  // SRC_IGA_BASIS_FUNCTION_HANDLER_H_
