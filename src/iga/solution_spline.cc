/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <vector>

#include "solution_spline.h"

iga::SolutionSpline::SolutionSpline(const std::shared_ptr<spl::NURBS<2>> &spl, const arma::dvec &solution) {
  std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector = {spl->GetKnotVector(0), spl->GetKnotVector(1)};
  std::array<Degree, 2> degree = {spl->GetDegree(0), spl->GetDegree(1)};
  std::vector<double> weights = spl->GetWeights();
  std::vector<baf::ControlPoint> control_points;
  std::vector<double> cp = spl->GetControlPoints();
  uint64_t l = 0;
  for (uint64_t i = 0; i < cp.size() - 2; i += 3) {
    control_points.emplace_back(baf::ControlPoint{cp[i], cp[i + 1], solution(l)});
    l += 1;
  }
  solution_spl_ = std::make_shared<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
}

std::shared_ptr<spl::NURBS<2>> iga::SolutionSpline::GetSolutionSpline() {
  return solution_spl_;
}
