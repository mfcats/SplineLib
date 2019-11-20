/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_BAF_B_SPLINE_BASIS_FUNCTION_H_
#define SRC_BAF_B_SPLINE_BASIS_FUNCTION_H_

#include "src/baf/knot_vector.h"
#include "src/util/named_type.h"

namespace splinelib::src::baf {
class BSplineBasisBasisFunction {
 public:
  BSplineBasisBasisFunction(BSplineBasisBasisFunction const &other) = delete;
  BSplineBasisBasisFunction(BSplineBasisBasisFunction &&other) = delete;
  virtual BSplineBasisBasisFunction & operator=(BSplineBasisBasisFunction const &rhs) = delete;
  virtual BSplineBasisBasisFunction & operator=(BSplineBasisBasisFunction &&rhs) = delete;
  virtual ~BSplineBasisBasisFunction() = default;

  // The evaluation of the i-th basis function of degree p > 0 N_{i,p} is a linear combination of the basis functions
  // N_{i,p-1} and N_{i+1,p-1} (see NURBS book equation 2.5). Therefore, for each basis function of degree > 0 a pointer
  // to these two basis functions is set in constructor, so that a basis function can be evaluated recursively.
  double Evaluate(const ParametricCoordinate &ParametricCoordinate) const;

  double EvaluateDerivative(const ParametricCoordinate &param_coord, const Derivative &derivative) const;

 protected:
  BSplineBasisBasisFunction(const KnotVector &knot_vector, const Degree &degree, const KnotSpan &start_of_support);

  virtual double EvaluateOnSupport(const ParametricCoordinate &param_coord) const = 0;
  virtual double EvaluateDerivativeOnSupport(const ParametricCoordinate &param_coord,
                                             const Derivative &derivative) const = 0;

  Degree GetDegree() const;
  ParametricCoordinate GetStartKnot() const;
  ParametricCoordinate GetEndKnot() const;

 private:
  bool IsCoordinateInSupport(const ParametricCoordinate &param_coord) const;

  Degree degree_;
  ParametricCoordinate start_knot_{0};
  ParametricCoordinate end_knot_{1};
  bool end_knot_is_last_knot_{false};
};
}  // namespace splinelib::src::baf

#endif  // SRC_BAF_B_SPLINE_BASIS_FUNCTION_H_
