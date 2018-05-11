/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_BASIS_FUNCTION_H_
#define SRC_BASIS_FUNCTION_H_

#include <memory>
#include <utility>
#include <vector>

#include "knot_vector.h"

namespace baf {
class BasisFunction {
 public:
  double Evaluate(double paramCoord) const;

  double EvaluateDerivative(int derivative, double param_coord) const;

 protected:
  BasisFunction(const KnotVector &knot_vector, int degree, uint64_t start);

  double GetKnot(uint64_t knot_position) const;

  uint64_t GetStartOfSupport() const;

  int GetDegree() const;

  virtual double EvaluateOnSupport(double param_coord) const = 0;

  virtual double EvaluateDerivativeOnSupport(int derivative, double param_coord) const = 0;

 private:
  bool IsCoordinateInSupport(double param_coord) const;

  bool IsCoordinateInSupportSpan(double param_coord) const;

  bool IsCoordinateSpecialCaseWithLastKnot(double param_coord) const;

  KnotVector knotVector_;
  int degree_;
  uint64_t start_of_support_;
};
}

#endif  // SRC_BASIS_FUNCTION_H_
