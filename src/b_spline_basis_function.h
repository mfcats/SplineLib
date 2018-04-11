/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_BSPLINEBASISFUNCTION_H
#define SPLINELIB_BSPLINEBASISFUNCTION_H

#include <cmath>
#include <limits>
#include <memory>

#include <iostream>

#include "basis_function.h"
#include "numeric_settings.h"

class BSplineBasisFunction : public BasisFunction {
 public:
  BSplineBasisFunction(KnotVector knot_vector,
                       Degree deg,
                       uint64_t start_of_support);

 protected:
  double EvaluateOnSupport(double param_coord) const override;
  double EvaluateDerivativeOnSupport(Derivative derivative,
                                     double param_coord) const override;

  std::unique_ptr<BasisFunction> left_lower_degree_;
  std::unique_ptr<BasisFunction> right_lower_degree_;

 private:
  void SetLowerDegreeBsisFunctions(KnotVector knot_vector,
                                   uint64_t start_of_support,
                                   Degree deg);

  double ComputeLeftQuotientDenominatorInverse() const;
  double ComputeRightQuotientDenominatorInverse() const;

  double ComputeLeftQuotient(double param_coord) const;
  double ComputeRightQuotient(double param_coord) const;
};

#endif // SPLINELIB_BSPLINEBASISFUNCTION_H
