/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_ELM_ELEMENT_INTEGRATION_POINT_H_
#define SRC_IGA_ELM_ELEMENT_INTEGRATION_POINT_H_

#include <array>
#include <vector>

namespace iga {
namespace elm {
class ElementIntegrationPoint {
 public:
  ElementIntegrationPoint(std::vector<double> basis_functions, double weight, double jac_det);

  ElementIntegrationPoint(std::array<std::vector<double>, 2> basis_function_derivatives, double weight, double jac_det);

  std::vector<double> GetNonZeroBasisFunctions() const;
  std::vector<double> GetNonZeroBasisFunctionDerivatives(int dir) const;
  double GetWeight() const;
  double GetJacobianDeterminant() const;
  int GetNumberOfNonZeroBasisFunctions() const;
  int GetNumberOfNonZeroBasisFunctionDerivatives(int dir) const;
  double GetBasisFunctionValue(int firstNonZeroOffset) const;
  double GetBasisFunctionDerivativeValue(int firstNonZeroOffset, int dir) const;

 private:
  std::vector<double> non_zero_basis_functions_;
  std::array<std::vector<double>, 2> non_zero_basis_function_derivatives_;
  double weight_;
  double jac_det_;
};
}  // namespace elm
}  // namespace iga

#endif  // SRC_IGA_ELM_ELEMENT_INTEGRATION_POINT_H_
