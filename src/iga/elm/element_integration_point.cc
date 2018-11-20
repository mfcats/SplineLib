/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "element_integration_point.h"

iga::elm::ElementIntegrationPoint::ElementIntegrationPoint(std::vector<double> basis_functions, double weight,
    double jac_det) : non_zero_basis_functions_(std::move(basis_functions)), weight_(weight), jac_det_(jac_det) {}

iga::elm::ElementIntegrationPoint::ElementIntegrationPoint(
    std::array<std::vector<double>, 2> basis_function_derivatives, double weight, double jac_det) :
    non_zero_basis_function_derivatives_(std::move(basis_function_derivatives)), weight_(weight), jac_det_(jac_det) {
}

std::vector<double> iga::elm::ElementIntegrationPoint::GetNonZeroBasisFunctions() const {
  return non_zero_basis_functions_;
}

std::vector<double> iga::elm::ElementIntegrationPoint::GetNonZeroBasisFunctionDerivatives(int dir) const {
  return non_zero_basis_function_derivatives_[dir];
}

double iga::elm::ElementIntegrationPoint::GetWeight() const {
  return weight_;
}

double iga::elm::ElementIntegrationPoint::GetJacobianDeterminant() const {
  return jac_det_;
}

int iga::elm::ElementIntegrationPoint::GetNumberOfNonZeroBasisFunctions() const {
  return static_cast<int>(non_zero_basis_functions_.size());
}

int iga::elm::ElementIntegrationPoint::GetNumberOfNonZeroBasisFunctionDerivatives(int dir) const {
  return static_cast<int>(non_zero_basis_function_derivatives_[dir].size());
}

double iga::elm::ElementIntegrationPoint::GetBasisFunctionValue(int firstNonZeroOffset) const {
#ifdef DEBUG
  return non_zero_basis_functions_.at(firstNonZeroOffset);
#else
  return non_zero_basis_functions_[firstNonZeroOffset];
#endif
}

double iga::elm::ElementIntegrationPoint::GetBasisFunctionDerivativeValue(int firstNonZeroOffset, int dir) const {
#ifdef DEBUG
  return non_zero_basis_function_derivatives_.at(dir).at(firstNonZeroOffset);
#else
  return non_zero_basis_function_derivatives_[dir][firstNonZeroOffset];
#endif
}
