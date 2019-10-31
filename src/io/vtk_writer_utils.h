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

#include "src/util/any_casts.h"
#include "src/spl/spline.h"

namespace splinelib::src::io {
template<int PARAMETRIC_DIMENSIONALITY>
class VTKWriterUtils {
 public:
  static std::array<double, 2 * PARAMETRIC_DIMENSIONALITY> GetEdgeKnots(
      std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>> spline_ptr) {
    std::array<double, 2 * PARAMETRIC_DIMENSIONALITY> knots{};
    for (int i = 0; i < 2 * PARAMETRIC_DIMENSIONALITY; ++i) {
      knots[i] = i < PARAMETRIC_DIMENSIONALITY ? spline_ptr->GetKnotVector(i)->GetKnot(0).Get()
                                               : spline_ptr->GetKnotVector(
              i - PARAMETRIC_DIMENSIONALITY)->GetLastKnot().Get();
    }
    return knots;
  }

  static int NumberOfCells(std::array<int, PARAMETRIC_DIMENSIONALITY> scattering) {
    int product = 1;
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      product *= scattering[i];
    }
    return product;
  }

  static std::array<int, PARAMETRIC_DIMENSIONALITY> GetPointHandlerLength(
      std::array<int, PARAMETRIC_DIMENSIONALITY> scattering) {
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      ++scattering[i];
    }
    return scattering;
  }

  static void WritePoints(std::ofstream &file,
                          const std::any &spline,
                          std::array<int, PARAMETRIC_DIMENSIONALITY> scattering) {
    std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>>
        spline_ptr = util::any_casts::GetSpline<PARAMETRIC_DIMENSIONALITY>(spline);
    std::array<double, 2 * PARAMETRIC_DIMENSIONALITY> knots = GetEdgeKnots(spline_ptr);
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(GetPointHandlerLength(scattering));
    for (int i = 0; i < point_handler.GetNumberOfTotalMultiIndices(); ++point_handler, ++i) {
      std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> coords{};
      for (int j = 0; j < PARAMETRIC_DIMENSIONALITY; ++j) {
        coords[j] = ParametricCoordinate(
            knots[j] + point_handler[Dimension{j}] * (knots[j + PARAMETRIC_DIMENSIONALITY] - knots[j]) / scattering[j]);
      }
      for (int k = 0; k < 3; ++k) {
        file << (k < spline_ptr->GetPointDim() ? spline_ptr->Evaluate(coords, {k})[0] : 0) << (k < 2 ? ' ' : '\n');
      }
    }
  }

  static void WriteCellTypes(std::ofstream &file,
                             std::array<int, PARAMETRIC_DIMENSIONALITY> scattering,
                             int cell_type) {
    for (int i = 0; i < NumberOfCells(scattering); ++i) {
      file << cell_type << "\n";
    }
  }
};
}  // namespace splinelib::src::io

#endif  // SRC_IO_VTK_WRITER_UTILS_H_
