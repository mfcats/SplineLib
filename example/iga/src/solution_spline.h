/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_SOLUTION_SPLINE_H_
#define SRC_IGA_SOLUTION_SPLINE_H_

#include <armadillo>
#include <vector>

#include "nurbs.h"

namespace iga {
template<int DIM>
class SolutionSpline {
 public:
  SolutionSpline(const std::shared_ptr<spl::NURBS<DIM>> &spl, const arma::dvec &solution) {
    std::vector<baf::ControlPoint> control_points;
    std::array<int, DIM> points_per_direction = spl->GetPointsPerDirection();
    util::MultiIndexHandler<DIM> point_handler(points_per_direction);
    for (uint64_t i = 0, l = 0; i < static_cast<uint64_t>(spl->GetNumberOfControlPoints()); ++i, ++point_handler, ++l) {
      std::vector<double> temp;
      for (int j = 0; j < spl->GetPointDim(); ++j) {
        temp.emplace_back(spl->GetControlPoint(point_handler.GetIndices(), j));
      }
      temp.emplace_back(solution(l));
      control_points.emplace_back(baf::ControlPoint(temp));
    }
    solution_spl_ = std::make_shared<spl::NURBS<DIM>>(*spl.get(), control_points);
  }

  std::shared_ptr<spl::NURBS<DIM>> GetSolutionSpline() {
    return solution_spl_;
  }

 private:
  std::shared_ptr<spl::NURBS<DIM>> solution_spl_;
};
}  // namespace iga

#endif  // SRC_IGA_SOLUTION_SPLINE_H_
