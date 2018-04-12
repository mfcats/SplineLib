/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_ZERODEGREEBSPLINEBASISFUNCTION_H
#define SPLINELIB_ZERODEGREEBSPLINEBASISFUNCTION_H

#include <vector>

#include "basis_function.h"
#include "knot_vector.h"

class ZeroDegreeBSplineBasisFunction : public BasisFunction {
 public:
  explicit ZeroDegreeBSplineBasisFunction(const KnotVector &knot_vector,
                                          uint64_t start_of_support);

 protected:
  double EvaluateOnSupport(double param_coord) const override;
  double EvaluateDerivativeOnSupport(Derivative derivative,
                                     double param_coord) const override;
};

#endif // SPLINELIB_ZERODEGREEBSPLINEBASISFUNCTION_H
