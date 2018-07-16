/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_ELM_ELEMENT_INTEGRATION_POINT_H_
#define SRC_ELM_ELEMENT_INTEGRATION_POINT_H_

#include <vector>

namespace elm {
class ElementIntegrationPoint {
 public:
  explicit ElementIntegrationPoint(const std::vector<double> &non_zero_basis_functions);

  std::vector<double> GetNonZeroBasisFunctions() const;
  int GetNumberOfNonZeroBasisFunctions() const;
  double GetBasisFunctionValue(int firstNonZeroOffset) const;

 private:
  std::vector<double> non_zero_basis_functions_;
};
}  // namespace elm

#endif  // SRC_ELM_ELEMENT_INTEGRATION_POINT_H_
