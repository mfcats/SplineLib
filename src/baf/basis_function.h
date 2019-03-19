/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_BAF_BASIS_FUNCTION_H_
#define SRC_BAF_BASIS_FUNCTION_H_

#include "knot_vector.h"
#include "named_type.h"

using Degree = util::NamedType<int, struct DegreeParameter>;
using Derivative = util::NamedType<int, struct DerivativeParameter>;

namespace baf {
class BasisFunction {
 public:
  virtual ~BasisFunction() = default;

  // The evaluation of the i-th basis function of degree p > 0 N_{i,p} is a linear combination of the basis functions
  // N_{i,p-1} and N_{i+1,p-1} (see NURBS book equation 2.5). Therefore, for each basis function of degree > 0 a pointer
  // to these two basis functions is set in constructor, so that a basis function can be evaluated recursively.
  double Evaluate(const ParamCoord &paramCoord) const;

  double EvaluateDerivative(const ParamCoord &param_coord, const Derivative &derivative) const;

 protected:
  BasisFunction(const KnotVector &knot_vector, const Degree &degree, const KnotSpan &start_of_support);

  virtual double EvaluateOnSupport(const ParamCoord &param_coord) const = 0;
  virtual double EvaluateDerivativeOnSupport(const ParamCoord &param_coord, const Derivative &derivative) const = 0;

  Degree GetDegree() const;
  ParamCoord GetStartKnot() const;
  ParamCoord GetEndKnot() const;

 private:
  bool IsCoordinateInSupport(const ParamCoord &param_coord) const;

  Degree degree_;
  ParamCoord start_knot_{0};
  ParamCoord end_knot_{1};
  bool end_knot_is_last_knot_{false};
};
}  // namespace baf

#endif  // SRC_BAF_BASIS_FUNCTION_H_
