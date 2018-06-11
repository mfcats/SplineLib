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

#include <memory>
#include <utility>
#include <vector>

#include "knot_vector.h"

namespace baf {
class BasisFunction {
 public:
  // The evaluation of the i-th basis function of degree p > 0 N_{i,p} is a linear combination of the basis functions
  // N_{i,p-1} and N_{i+1,p-1} (see NURBS book equation 2.5). Therefore, for each basis function of degree > 0 a pointer
  // to these two basis functions is set in constructor, so that a basis function can be evaluated recursively.
  double Evaluate(ParamCoord paramCoord) const;

  double EvaluateDerivative(ParamCoord param_coord, int derivative) const;

 protected:
  BasisFunction(const KnotVector &knot_vector, int degree, uint64_t start);

  ParamCoord GetKnot(uint64_t knot_position) const;

  uint64_t GetStartOfSupport() const;

  int GetDegree() const;

  virtual double EvaluateOnSupport(ParamCoord param_coord) const = 0;

  virtual double EvaluateDerivativeOnSupport(ParamCoord param_coord, int derivative) const = 0;

 private:
  // Check if parametric coordinate is in knot vector range (see IsCoordinateInSupportSpan) and
  // if it is either in support span or meets the special case (see IsCoordinateSpecialCaseWithLastKnot).
  bool IsCoordinateInSupport(ParamCoord param_coord) const;

  // Check if parametric coordinate is in the range of knot spans where the basis function is defined to be non-zero.
  bool IsCoordinateInSupportSpan(ParamCoord param_coord) const;

  // Check if parametric coordinate is last knot of knot vector and last knot of basis function support range.
  bool IsCoordinateSpecialCaseWithLastKnot(ParamCoord param_coord) const;

  KnotVector knotVector_;
  int degree_;
  uint64_t start_of_support_;
};
}

#endif  // SRC_BAF_BASIS_FUNCTION_H_
