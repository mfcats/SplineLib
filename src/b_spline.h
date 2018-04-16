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

#include <vector>

#include "control_point.h"
#include "parameter_space.h"

class BSpline {
 public:
  BSpline(const KnotVector &knot_vector,
          Degree degree,
          const std::vector<ControlPoint> &control_points);

  std::vector<double> Evaluate(double param_coord, std::vector<Dimension> dimensions) const;
  std::vector<double> EvaluateDerivative(double param_coord,
                                         std::vector<int> dimensions,
                                         int derivative) const;

  int degree() const;
  KnotVector knot_vector() const;

 private:
  std::vector<double> ExtractControlPointValues(double param_coord, int dimension) const;
  double ComputeWeightedSum(std::vector<double> basis_function_values,
                            std::vector<double> control_point_values) const;

  ParameterSpace parameter_space_;
  std::vector<std::vector<double>> control_points_;
};

#endif // SPLINELIB_B_SPLINE_H
