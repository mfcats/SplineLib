/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_PARAMETER_SPACE_H_
#define SRC_SPL_PARAMETER_SPACE_H_

#include <vector>

#include "basis_function.h"
#include "element.h"
#include "integration_rule.h"
#include "knot_vector.h"
#include "element_integration_point.h"

namespace spl {
class ParameterSpace {
 public:
  ParameterSpace() {}

  ParameterSpace(const baf::KnotVector &knot_vector, int degree);

  std::vector<double> EvaluateAllNonZeroBasisFunctions(ParamCoord param_coord) const;
  std::vector<double> EvaluateAllNonZeroBasisFunctionDerivatives(ParamCoord param_coord, int derivative) const;

  std::vector<std::unique_ptr<baf::BasisFunction>>::const_iterator GetFirstNonZeroKnot(ParamCoord param_coord) const;
  int degree() const;

  baf::KnotVector knot_vector() const;

  std::vector<elm::ElementIntegrationPoint>
  EvaluateAllElementNonZeroBasisFunctions(int element_number, const itg::IntegrationRule<1> &rule) const;
  std::vector<elm::ElementIntegrationPoint>
  EvaluateAllElementNonZeroBasisFunctionDerivatives(int element_number, const itg::IntegrationRule<1> &rule) const;

  std::vector<elm::Element> GetElementList() const;
  ParamCoord ReferenceSpace2ParameterSpace(double upper, double lower, double point) const;

 private:
  baf::KnotVector knot_vector_;
  int degree_;
  std::vector<std::unique_ptr<baf::BasisFunction>> basis_functions_;
};
}  // namespace spl

#endif  // SRC_SPL_PARAMETER_SPACE_H_
