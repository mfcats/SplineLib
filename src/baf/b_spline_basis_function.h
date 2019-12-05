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

// BSplineBasisFunctions are piecewise polynomial functions of a specific degree forming a basis for the vector space of
// all piecewise polynomial functions of the desired degree for a fixed knot vector sequence. Each BSplineBasisFunction
// N_{i,p} is only nonzero on a limited number of subintervals of the knot vector.
// Therefore, not the entire knot vector is saved but only the interval where the BSplineBasisFunction is nonzero. As
// computation for the last knot is different, it is also stored if the end knot is the last knot of the knot vector.
namespace splinelib::src::baf {
class BSplineBasisFunction {
 public:
  using ConstBSplineBasisFunctionIterator = std::vector<std::shared_ptr<BSplineBasisFunction>>::const_iterator;

  // Computation of a basis functions N_{i,p} requires specification of a knot vector, a knot span i where local support
  // starts and a degree p. This function also creates all basis functions of lower degree required for evaluation.
  static BSplineBasisFunction * CreateDynamic(KnotVector const &knot_vector, KnotSpan const &start_of_support,
                                              Degree const &degree);

  BSplineBasisFunction() = default;
  BSplineBasisFunction(BSplineBasisFunction const &other) = delete;
  BSplineBasisFunction(BSplineBasisFunction &&other) = delete;
  virtual BSplineBasisFunction & operator=(BSplineBasisFunction const &rhs) = delete;
  BSplineBasisFunction & operator=(BSplineBasisFunction &&rhs) = delete;
  virtual ~BSplineBasisFunction() = default;

  // The recurrence evaluation formula due to DeBoor, Cox, and Mansfield is implemented (see NURBS book equation 2.5).
  double Evaluate(ParametricCoordinate const &parametric_coordinate) const;
  // For derivative evaluation the formula 2.9 in NURBS book is implemented. For 0th derivative Evaluate is called and
  // for parametric coordinate not in local support zero is returned.
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
