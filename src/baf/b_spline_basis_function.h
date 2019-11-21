/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_BAF_B_SPLINE_BASIS_FUNCTION_H_
#define SRC_BAF_B_SPLINE_BASIS_FUNCTION_H_

#include <vector>

#include "src/baf/knot_vector.h"
#include "src/util/named_type.h"

namespace splinelib::src::baf {
class BSplineBasisFunction {
 public:
  using ConstBSplineBasisFunctionIterator = std::vector<std::shared_ptr<BSplineBasisFunction>>::const_iterator;

  static BSplineBasisFunction * CreateDynamic(KnotVector const &knot_vector, KnotSpan const &start_of_support,
                                              Degree const &degree);

  BSplineBasisFunction() = default;
  BSplineBasisFunction(BSplineBasisFunction const &other) = delete;
  BSplineBasisFunction(BSplineBasisFunction &&other) noexcept;
  virtual BSplineBasisFunction & operator=(BSplineBasisFunction const &rhs) = delete;
  BSplineBasisFunction & operator=(BSplineBasisFunction &&rhs) noexcept;
  virtual ~BSplineBasisFunction() = default;

  // The evaluation of the i-th basis function of degree p > 0 N_{i,p} is a linear combination of the basis functions
  // N_{i,p-1} and N_{i+1,p-1} (see NURBS book equation 2.5). Therefore, for each basis function of degree > 0 a pointer
  // to these two basis functions is set in constructor, so that a basis function can be evaluated recursively.
  double Evaluate(ParametricCoordinate const &parametric_coordinate) const;
  double EvaluateDerivative(ParametricCoordinate const &parametric_coordinate, Derivative const &derivative) const;

 protected:
  BSplineBasisFunction(KnotVector const &knot_vector, KnotSpan const &start_of_support, Degree const &degree);

  Degree GetDegree() const;
  ParametricCoordinate GetStartKnot() const;
  ParametricCoordinate GetEndKnot() const;

  virtual double EvaluateOnSupport(ParametricCoordinate const &parametric_coordinate) const = 0;
  virtual double EvaluateDerivativeOnSupport(ParametricCoordinate const &parametric_coordinate,
                                             Derivative const &derivative) const = 0;

 private:
  bool IsCoordinateInSupport(ParametricCoordinate const &parametric_coordinate) const;

  Degree degree_{};
  ParametricCoordinate start_knot_{0};
  ParametricCoordinate end_knot_{1};
  bool end_knot_is_last_knot_{false};
};

#include "src/baf/b_spline_basis_funtion.inc"
}  // namespace splinelib::src::baf

#endif  // SRC_BAF_B_SPLINE_BASIS_FUNCTION_H_
