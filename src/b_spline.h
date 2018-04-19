/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_B_SPLINE_H_
#define SRC_B_SPLINE_H_

#include <vector>

#include "control_point.h"
#include "parameter_space.h"
#include "element_generator.h"
#include "knot_vector.h"

class BSpline {
 public:
  BSpline(const KnotVector &knot_vector, int degree, const std::vector<ControlPoint> &control_points);

  std::vector<double> Evaluate(double param_coord, const std::vector<int> &dimensions) const;
  std::vector<double> EvaluateDerivative(double param_coord, const std::vector<int> &dimensions, int derivative) const;

  int GetDegree() const;
  KnotVector GetKnotVector() const;

  std::vector<Element> GetElementList();

  std::vector<std::unique_ptr<BasisFunction>> *basis_functions() {
    return &basis_functions_;
  }

 private:
  std::vector<double> ExtractControlPointValues(double param_coord, int dimension) const;
  double ComputeWeightedSum(const std::vector<double> &basis_function_values,
                            std::vector<double> control_point_values) const;

  ParameterSpace parameter_space_;
  std::vector<double> control_points_;
  int dim;
};

#endif  // SRC_B_SPLINE_H_
