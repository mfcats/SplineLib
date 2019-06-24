/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University
This file is part of SplineLib.
SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.
SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_SOLUTION_VTK_WRITER_H_
#define SRC_IGA_SOLUTION_VTK_WRITER_H_

#include <armadillo>
#include <fstream>
#include <string>
#include <vector>

#include "multi_index_handler.h"
#include "nurbs.h"
#include "solution_spline.h"
#include "vtk_writer.h"

namespace iga {
template<int DIM>
class SolutionVTKWriter {
 public:
  SolutionVTKWriter() = default;

  void WriteSolutionToVTK(const std::shared_ptr<spl::NURBS<DIM>> &spl, const arma::dvec &solution,
                          const std::vector<std::vector<int>> &scattering, const std::string &filename) const {
    iga::SolutionSpline<DIM> sol_spl(spl, solution);
    std::shared_ptr<spl::NURBS<DIM>> solution_spl = sol_spl.GetSolutionSpline();
    int cp_dim = spl->GetPointDim();
    std::vector<double> dxi;
    std::array<int, DIM> num_pnts{};
    for (int i = 0; i < DIM; ++i) {
      baf::KnotVector knots = *spl->GetKnotVector(i);
      dxi.emplace_back((knots.GetKnot(knots.GetNumberOfKnots() - 1) - knots.GetKnot(0)).get() / scattering[0][i]);
      num_pnts[i] = scattering[0][i] + 1;
    }
    std::vector<double> point_data;
    util::MultiIndexHandler<DIM> mih(num_pnts);
    while (true) {
      std::array<ParamCoord, DIM> param_coords{};
      for (int i = 0; i < DIM; ++i) {
        param_coords[i] = spl->GetKnotVector(i)->GetKnot(0) + ParamCoord{mih[i] * dxi[i]};
      }
      point_data.emplace_back(solution_spl->Evaluate(param_coords, {cp_dim})[0]);
      if (mih.Get1DIndex() == mih.Get1DLength() - 1) break;
      ++mih;
    }
    io::VTKWriter vtk_writer;
    std::vector<std::any> splines = {std::make_any<std::shared_ptr<spl::NURBS<DIM>>>(spl)};
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
};
}  // namespace iga

#endif  // SRC_IGA_SOLUTION_VTK_WRITER_H_
