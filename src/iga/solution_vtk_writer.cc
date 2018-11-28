/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "solution_vtk_writer.h"

iga::SolutionVTKWriter::SolutionVTKWriter() = default;

void iga::SolutionVTKWriter::WriteSolutionToVTK(const std::shared_ptr<spl::NURBS<2>> &spl, const arma::dvec &solution,
    const std::vector<std::vector<int>> &scattering, const std::string &filename) {
  KnotVectors<2> knot_vector = {spl->GetKnotVector(0), spl->GetKnotVector(1)};
  std::array<Degree, 2> degree = {spl->GetDegree(0), spl->GetDegree(1)};
  std::vector<double> weights = spl->GetWeights();
  std::vector<baf::ControlPoint> control_points;
  std::vector<double> cp = spl->GetControlPoints();
  uint64_t l = 0;
  for (uint64_t i = 0; i < cp.size() - 2; i += 3) {
    control_points.emplace_back(baf::ControlPoint{cp[i], cp[i + 1], solution(l)});
    l += 1;
  }
  spl::NURBS<2> solution_spl(knot_vector, degree, control_points, weights);
  double dxi = (spl->GetKnots()[0][spl->GetKnots()[0].size() - 1] - spl->GetKnots()[0][0]).get() / scattering[0][0];
  double deta = (spl->GetKnots()[1][spl->GetKnots()[1].size() - 1] - spl->GetKnots()[1][0]).get() / scattering[0][1];
  std::vector<double> point_data;
  std::array<ParamCoord, 2> param_coords = {spl->GetKnots()[0][0], spl->GetKnots()[1][0]};
  for (int j = 0; j <= scattering[0][1]; ++j) {
    for (int i = 0; i <= scattering[0][0]; ++i) {
      point_data.emplace_back(solution_spl.Evaluate(param_coords, {2})[0]);
      param_coords[0] = param_coords[0] + ParamCoord{dxi};
    }
    param_coords[1] = param_coords[1] + ParamCoord{deta};
    param_coords[0] = spl->GetKnots()[0][0];
  }
  io::VTKWriter vtk_writer;
  std::vector<std::any> splines = {std::make_any<std::shared_ptr<spl::NURBS<2>>>(spl)};
  vtk_writer.WriteFile(splines, filename, scattering);
  std::ofstream newFile;
  newFile.open(filename, std::ofstream::out | std::ofstream::app);
  if (newFile.is_open()) {
    newFile << "\nPOINT_DATA " << point_data.size() << "\nSCALARS solution float 1\nLOOKUP_TABLE default\n";
    for (auto &p : point_data) {
      newFile << p << "\n";
    }
    newFile.close();
  }
}
