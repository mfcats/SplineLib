/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_B_SPLINE_H
#define SPLINELIB_B_SPLINE_H

#include <utility>
#include <vector>

#include "basis_function.h"
#include "control_point.h"
#include "knot_vector.h"

class BSpline {
 public:
  BSpline(KnotVector knot_vector,
          Degree degree,
          std::vector<ControlPoint> control_points);

  std::vector<double> Evaluate(double param_coord, std::vector<Dimension> dimensions) const;
  std::vector<double> EvaluateDerivative(double param_coord,
                                         std::vector<int> dimensions,
                                         int derivative) const;

  int degree() const {
    return degree_;
  }
  KnotVector knot_vector() const {
    return knot_vector_;
  }
  std::vector<std::unique_ptr<BasisFunction>> *basis_functions() {
    return &basis_functions_;
  }

 private:
  std::vector<double> EvaluateAllNonZeroBasisFunctions(double param_coord) const;
  std::vector<double> EvaluateAllNonZeroBasisFunctionDerivatives(double param_coord,
                                                                 int derivative) const;
  std::vector<double> ExtractControlPointValues(double param_coord, int dimension) const;
  double ComputeWeightedSum(std::vector<double> basis_function_values,
                            std::vector<double> control_point_values) const;

  KnotVector knot_vector_;
  Degree degree_;
  std::vector<ControlPoint> control_points_;
  std::vector<std::unique_ptr<BasisFunction>> basis_functions_;
};

#endif // SPLINELIB_B_SPLINE_H
