/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_PARAMETER_SPACE_H_
#define SRC_PARAMETER_SPACE_H_

#include <vector>

#include "basis_function.h"
#include "knot_vector.h"
#include "integration_rule.h"
#include "element.h"

class ParameterSpace {
 public:
  ParameterSpace(const KnotVector &knot_vector, int degree);

  std::vector<double> EvaluateAllNonZeroBasisFunctions(double param_coord) const;
  std::vector<double> EvaluateAllNonZeroBasisFunctionDerivatives(double param_coord, int derivative) const;
  int degree() const;
  KnotVector knot_vector() const;
  std::vector<std::vector<double>> EvaluateAllElementNonZeroBasisFunctions(int element_number, IntegrationRule<1> rule);
  std::vector<Element> GetElementList();

 private:
  std::vector<std::unique_ptr<BasisFunction>>::const_iterator GetFirstNonZeroKnot(double param_coord) const;
  double TransformElementPoint(double upper, double lower, double point) const;

  KnotVector knot_vector_;
  int degree_;
  std::vector<std::unique_ptr<BasisFunction>> basis_functions_;
};

#endif  // SRC_PARAMETER_SPACE_H_
