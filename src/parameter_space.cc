/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "parameter_space.h"

#include "basis_function_factory.h"

ParameterSpace::ParameterSpace(const KnotVector &knot_vector, int degree) : degree_(degree), knot_vector_(knot_vector) {
  BasisFunctionFactory factory;

  basis_functions_.reserve(knot_vector_.Size() - degree_ - 1);
  for (uint64_t i = 0; i < (knot_vector_.Size() - degree_ - 1); ++i) {
    basis_functions_.emplace_back(factory.CreateDynamic(knot_vector_, i, degree_));
  }
}

std::vector<double> ParameterSpace::EvaluateAllNonZeroBasisFunctions(double param_coord) const {
  auto first_non_zero = GetFirstNonZeroKnot(param_coord);
  std::vector<double> basis_function_values(static_cast<u_int64_t >(degree_) + 1, 0.0);
  for (int i = 0; i < degree_ + 1; ++i) {
    basis_function_values[i] = (*first_non_zero)->Evaluate(param_coord);
    ++first_non_zero;
  }
  return basis_function_values;
}

std::vector<double> ParameterSpace::EvaluateAllNonZeroBasisFunctionDerivatives(double param_coord,
                                                                               int derivative) const {
  auto first_non_zero = GetFirstNonZeroKnot(param_coord);
  std::vector<double> basis_function_values(static_cast<u_int64_t >(degree_) + 1, 0.0);
  for (int i = 0; i < degree_ + 1; ++i) {
    basis_function_values[i] = (*first_non_zero)->EvaluateDerivative(derivative, param_coord);
    ++first_non_zero;
  }
  return basis_function_values;
}

int ParameterSpace::degree() const {
  return degree_;
}

KnotVector ParameterSpace::knot_vector() const {
  return knot_vector_;
}

std::vector<std::unique_ptr<BasisFunction>>::const_iterator ParameterSpace::GetFirstNonZeroKnot(
    double param_coord) const {
  return basis_functions_.begin() + knot_vector_.GetKnotSpan(param_coord) - degree_;
}
