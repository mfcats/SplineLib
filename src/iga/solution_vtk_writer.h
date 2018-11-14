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

#include "nurbs.h"
#include "vtk_writer.h"

namespace iga {
class SolutionVTKWriter {
 public:
    SolutionVTKWriter();

    void WriteSolutionToVTK(const std::shared_ptr<spl::NURBS<2>> &spl, const arma::dvec &solution,
                            const std::vector<std::vector<int>> &scattering, const std::string &filename);
};
}  // namespace iga

#endif  // SRC_IGA_SOLUTION_VTK_WRITER_H_
