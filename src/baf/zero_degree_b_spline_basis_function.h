/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_BAF_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_H_
#define SRC_BAF_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_H_

#include "basis_function.h"
#include "knot_vector.h"

namespace baf {
class ZeroDegBSplBasFnc : public baf::BasisFunction {
 public:
  ZeroDegBSplBasFnc(const baf::KnotVector &knot_vector, const KnotSpan &start_of_support);

 protected:
  double EvaluateOnSupport(const ParamCoord &param_coord) const override;
  double EvaluateDerivativeOnSupport(const ParamCoord &param_coord, const Derivative &derivative) const override;
};
}  // namespace baf

#endif  // SRC_BAF_ZERO_DEGREE_B_SPLINE_BASIS_FUNCTION_H_
