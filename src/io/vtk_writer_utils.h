/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_VTK_WRITER_UTILS_H_
#define SRC_IO_VTK_WRITER_UTILS_H_

#include <array>
#include <fstream>

#include "any_casts.h"
#include "spline.h"

namespace splinelib::src::io {
template<int DIM>
class VTKWriterUtils {
 public:
  static std::array<double, 2 * DIM> GetEdgeKnots(std::shared_ptr<spl::Spline<DIM>> spline_ptr) {
    std::array<double, 2 * DIM> knots{};
    for (int i = 0; i < 2 * DIM; ++i) {
      knots[i] = i < DIM ? spline_ptr->GetKnotVector(i)->GetKnot(0).get()
                         : spline_ptr->GetKnotVector(i - DIM)->GetLastKnot().get();
    }
    return knots;
  }

  static int NumberOfCells(std::array<int, DIM> scattering) {
    int product = 1;
    for (int i = 0; i < DIM; ++i) {
      product *= scattering[i];
    }
    return product;
  }

  static std::array<int, DIM> GetPointHandlerLength(std::array<int, DIM> scattering) {
    for (int i = 0; i < DIM; ++i) {
      ++scattering[i];
    }
    return scattering;
  }

  static void WritePoints(std::ofstream &file, const std::any &spline, std::array<int, DIM> scattering) {
    std::shared_ptr<spl::Spline<DIM>> spline_ptr = util::AnyCasts::GetSpline<DIM>(spline);
    std::array<double, 2 * DIM> knots = GetEdgeKnots(spline_ptr);
    util::MultiIndexHandler<DIM> point_handler(GetPointHandlerLength(scattering));
    for (int i = 0; i < point_handler.Get1DLength(); ++point_handler, ++i) {
      std::array<baf::ParamCoord, DIM> coords{};
      for (int j = 0; j < DIM; ++j) {
        coords[j] = baf::ParamCoord(knots[j] + point_handler[j] * (knots[j + DIM] - knots[j]) / scattering[j]);
      }
      for (int k = 0; k < 3; ++k) {
        file << (k < spline_ptr->GetPointDim() ? spline_ptr->Evaluate(coords, {k})[0] : 0) << (k < 2 ? " " : "\n");
      }
    }
  }

  static void WriteCellTypes(std::ofstream &file, std::array<int, DIM> scattering, int cell_type) {
    for (int i = 0; i < NumberOfCells(scattering); ++i) {
      file << cell_type << "\n";
    }
  }
};
}  // namespace splinelib::src::splinelib::src::io

#endif  // SRC_IO_VTK_WRITER_UTILS_H_
