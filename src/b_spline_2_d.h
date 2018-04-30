/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_B_SPLINE_2_D_H
#define SPLINELIB_B_SPLINE_2_D_H

#include <vector>

#include "control_point.h"
#include "parameter_space.h"

class BSpline2D {
 public:
  BSpline2D(const std::vector<KnotVector> &knot_vectors,
            const std::vector<int> &degrees,
            const std::vector<std::vector<ControlPoint>> &control_points) : dim(control_points[0][0].GetDimension()) {
    for (int space = 0; space < knot_vectors.size(); space++) {
      parameter_spaces_.emplace_back(ParameterSpace(knot_vectors[space], degrees[space]));
    }
    for (auto &dimension : control_points) {
      for (auto &cp : dimension) {
        for (int i = 0; i < dim; ++i) {
          control_points_.emplace_back(std::vector<double>({cp.GetValue(i)}));
        }
      }
    }
  }

  std::vector<double> Evaluate(const std::vector<double> &param_coord, const std::vector<int> &dimensions) const {
    return {0};
  }
  std::vector<double> EvaluateDerivative(std::vector<double> param_coord,
                                         const std::vector<int> &dimensions,
                                         int derivative) const;

 private:
  std::vector<ParameterSpace> parameter_spaces_;
  std::vector<std::vector<double>> control_points_;
  int dim;
};

#endif //SPLINELIB_B_SPLINE_2_D_H
