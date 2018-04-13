/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_PARAMETER_SPACE_H
#define SPLINELIB_PARAMETER_SPACE_H

#include <vector>

#include "basis_function.h"
#include "knot_vector.h"

class ParameterSpace {
 public:
  ParameterSpace(const KnotVector &knot_vector, Degree degree);

  std::vector<double> EvaluateAllNonZeroBasisFunctions(double param_coord) const;
  std::vector<double> EvaluateAllNonZeroBasisFunctionDerivatives(double param_coord,
                                                                 int derivative) const;
  int degree() const;
  KnotVector knot_vector() const;
  std::vector<std::unique_ptr<BasisFunction>> *basis_functions();

 private:
  KnotVector knot_vector_;
  Degree degree_;
  std::vector<std::unique_ptr<BasisFunction>> basis_functions_;
};

#endif // SPLINELIB_B_SPLINE_H
