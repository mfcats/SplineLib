/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_VTK_WRITER_H_
#define SRC_IO_VTK_WRITER_H_

#include <any>
#include <fstream>
#include <string>
#include <vector>

#include "any_casts.h"

namespace io {
class VTKWriter {
 public:
  VTKWriter() = default;

  void WriteFile(const std::vector<std::any> &splines,
                 const std::string &filename,
                 const std::vector<std::vector<int>> &scattering) const {
    std::ofstream newFile;
    newFile.open(filename);
    if (newFile.is_open()) {
      newFile << "# vtk DataFile Version 3.0\nSpline from Splinelib\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";
      for (auto i = 0u; i < splines.size(); ++i) {
        AddSpline(newFile, splines[i], scattering[i]);
      }
      newFile.close();
    }
  }

 private:
  void AddSpline(std::ofstream &file, const std::any &spline, const std::vector<int> &scattering) const {
    int spline_dimension = util::AnyCasts::GetSplineDimension(spline);
    if (spline_dimension == 1) {
      Write1DSpline(file, spline, scattering[0]);
    } else if (spline_dimension == 2) {
      Write2DSpline(file, spline, {scattering[0], scattering[1]});
    } else if (spline_dimension == 3) {
      Write3DSpline(file, spline, {scattering[0], scattering[1], scattering[2]});
    } else {
      throw std::runtime_error("Only splines of dimensions 1 to 3 can be written to a vtk file.");
    }
  }

  void Write1DSpline(std::ofstream &file, const std::any &spline, int scattering) const {
    std::shared_ptr<spl::Spline<1>> spline_ptr = util::AnyCasts::GetSpline<1>(spline);
    WritePoints<1>(file, spline_ptr, {scattering});
    file << "\nCELLS " << scattering << " " << 3 * scattering << "\n";
    for (int i = 0; i < scattering; ++i) {
      file << "2 " << i << " " << i + 1 << "\n";
    }
    WriteCellTypes<1>(file, {scattering}, 3);
  }

  void Write2DSpline(std::ofstream &file, const std::any &spline, std::array<int, 2> scattering) const {
    std::shared_ptr<spl::Spline<2>> spline_ptr = util::AnyCasts::GetSpline<2>(spline);
    WritePoints<2>(file, spline_ptr, scattering);
    file << "\nCELLS " << NumberOfCells<2>(scattering) << " " << 5 * NumberOfCells<2>(scattering) << "\n";
    util::MultiIndexHandler<2> point_handler({scattering[0] + 1, scattering[1] + 1});
    for (int j = 0; j < scattering[0]; ++j) {
      for (int i = 0; i < scattering[1]; ++i) {
        file << "4 " << (scattering[0] + 1) * j + i << " " << (scattering[0] + 1) * j + i + 1 << " " <<
             (scattering[0] + 1) * (j + 1) + i + 1 << " " << (scattering[0] + 1) * (j + 1) + i << "\n";
      }
    }
    WriteCellTypes<2>(file, scattering, 9);
  }

  void Write3DSpline(std::ofstream &file, const std::any &spline, std::array<int, 3> scattering) const {
    std::shared_ptr<spl::Spline<3>> spline_ptr = util::AnyCasts::GetSpline<3>(spline);
    WritePoints<3>(file, spline_ptr, scattering);
    file << "\nCELLS " << NumberOfCells<3>(scattering) << " " << 9 * NumberOfCells<3>(scattering) << "\n";
    util::MultiIndexHandler<3> point_handler(GetPointHandlerLength<3>(scattering));
    for (; point_handler.Get1DIndex() < point_handler.Get1DLength() - 1;
           ++point_handler) {
      if (point_handler.GetIndices()[0] != scattering[0] && point_handler.GetIndices()[1] != scattering[1]
          && point_handler.GetIndices()[2] != scattering[2]) {
        file << "8 " << point_handler.Get1DIndex() << " " << (point_handler + 1).Get1DIndex() << " "
             << (point_handler + scattering[0] + 1).Get1DIndex() << " " << (point_handler - 1).Get1DIndex() << " "
             << (point_handler + scattering[1] * (scattering[0] + 1)).Get1DIndex() << " "
             << (point_handler + 1).Get1DIndex() << " " << (point_handler + scattering[0] + 1).Get1DIndex() << " "
             << (point_handler - 1).Get1DIndex() << std::endl;
        point_handler - (scattering[0] + 1) * (scattering[1] + 2);
      }
    }
    WriteCellTypes<3>(file, scattering, 12);
  }

  template<int DIM>
  int NumberOfCells(std::array<int, DIM> scattering) const {
    int product = 1;
    for (int i = 0; i < DIM; ++i) {
      product *= scattering[i];
    }
    return product;
  }

  template<int DIM>
  std::array<double, 2 * DIM> GetEdgeKnots(std::shared_ptr<spl::Spline<DIM>> spline_ptr) const {
    std::array<double, 2 * DIM> knots{};
    for (int i = 0; i < 2 * DIM; ++i) {
      knots[i] = i < DIM ? spline_ptr->GetKnotVector(i)->GetKnot(0).get()
                         : spline_ptr->GetKnotVector(i - DIM)->GetLastKnot().get();
    }
    return knots;
  }

  template<int DIM>
  std::array<int, DIM> GetPointHandlerLength(std::array<int, DIM> scattering) const {
    for (int i = 0; i < DIM; ++i) {
      ++scattering[i];
    }
    return scattering;
  }

  template<int DIM>
  void WritePoints(std::ofstream &file,
                   std::shared_ptr<spl::Spline<DIM>> spline_ptr,
                   std::array<int, DIM> scattering) const {
    file << "POINTS " << NumberOfCells<DIM>(GetPointHandlerLength<DIM>(scattering)) << " double\n";
    std::array<double, 2 * DIM> knots = GetEdgeKnots<DIM>(spline_ptr);
    util::MultiIndexHandler<DIM> point_handler(GetPointHandlerLength<DIM>(scattering));
    for (int i = 0; i < point_handler.Get1DLength(); ++point_handler, ++i) {
      std::array<ParamCoord, DIM> coords{};
      for (int j = 0; j < DIM; ++j) {
        coords[j] = ParamCoord(knots[j] + point_handler[j] * (knots[j + DIM] - knots[j]) / scattering[j]);
      }
      for (int k = 0; k < 3; ++k) {
        file << (k < spline_ptr->GetDimension() ? spline_ptr->Evaluate(coords, {k})[0] : 0) << (k < 2 ? " " : "\n");
      }
    }
  }

  template<int DIM>
  void WriteCellTypes(std::ofstream &file, std::array<int, DIM> scattering, int cell_type) const {
    file << "\nCELL_TYPES " << NumberOfCells<DIM>(scattering) << "\n";
    for (int j = 0; j < NumberOfCells<DIM>(scattering); ++j) {
      file << cell_type << "\n";
    }
  }
};
}  // namespace io
#endif  // SRC_IO_VTK_WRITER_H_
