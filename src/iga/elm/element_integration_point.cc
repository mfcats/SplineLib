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

#include <vector>

iga::elm::ElementIntegrationPoint::ElementIntegrationPoint(std::vector<double> basis_functions, double weight)
    : non_zero_basis_functions_(std::move(basis_functions)), weight_(weight), jac_det_(1.0) {}

iga::elm::ElementIntegrationPoint::ElementIntegrationPoint(std::vector<double> basis_functions, double weight,
    double jac_det) : non_zero_basis_functions_(std::move(basis_functions)), weight_(weight), jac_det_(jac_det) {}

std::vector<double> iga::elm::ElementIntegrationPoint::GetNonZeroBasisFunctions() const {
  return non_zero_basis_functions_;
}

double iga::elm::ElementIntegrationPoint::GetWeight() const {
  return weight_;
}

double iga::elm::ElementIntegrationPoint::GetJacDet() const {
  return jac_det_;
}

int iga::elm::ElementIntegrationPoint::GetNumberOfNonZeroBasisFunctions() const {
  return static_cast<int>(non_zero_basis_functions_.size());
}

double iga::elm::ElementIntegrationPoint::GetBasisFunctionValue(int firstNonZeroOffset) const {
#ifdef DEBUG
  return non_zero_basis_functions_.at(firstNonZeroOffset);
#else
  return non_zero_basis_functions_[firstNonZeroOffset];
#endif
}
