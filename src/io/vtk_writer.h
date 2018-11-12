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
#include "vtk_writer_utils.h"

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
      std::vector<int> dimensions = GetSplineDimensions(splines);
      std::vector<int> cells = GetNumberOfAllCells(dimensions, scattering);
      std::vector<int> points = GetNumberOfAllPoints(dimensions, scattering);
      newFile << "POINTS " << std::accumulate(points.begin(), points.end(), 0) << " double\n";
      for (auto i = 0u; i < splines.size(); ++i) {
        AddPoints(newFile, splines[i], scattering[i], dimensions[i]);
      }
      newFile << "\nCELLS " << std::accumulate(cells.begin(), cells.end(), 0) << " "
              << GetNumberOfCellEntries(dimensions, cells) << "\n";
      for (auto i = 0u; i < splines.size(); ++i) {
        AddCells(newFile, scattering[i], std::accumulate(points.begin(), points.begin() + i, 0), dimensions[i]);
      }
      newFile << "\nCELL_TYPES " << std::accumulate(cells.begin(), cells.end(), 0) << "\n";
      for (auto i = 0u; i < splines.size(); ++i) {
        AddCellTypes(newFile, scattering[i], dimensions[i]);
      }
      newFile.close();
    }
  }

 private:
  std::vector<int> GetSplineDimensions(const std::vector<std::any> &splines) const {
    std::vector<int> dimensions;
    for (const auto &spline : splines) {
      dimensions.push_back(util::AnyCasts::GetSplineDimension(spline));
    }
    return dimensions;
  }

  std::vector<int> GetNumberOfAllPoints(const std::vector<int> &dimensions,
                                        const std::vector<std::vector<int>> &scattering) const {
    std::vector<int> points;
    for (auto i = 0u; i < dimensions.size(); ++i) {
      points.push_back(dimensions[i] == 1 ? VTKWriterUtils<1>::NumberOfCells({scattering[i][0] + 1}) :
                       (dimensions[i] == 2 ? VTKWriterUtils<2>::NumberOfCells({scattering[i][0] + 1,
                                                                               scattering[i][1] + 1}) :
                        VTKWriterUtils<3>::NumberOfCells({scattering[i][0] + 1, scattering[i][1] + 1,
                                                          scattering[i][2] + 1})));
    }
    return points;
  }

  std::vector<int> GetNumberOfAllCells(const std::vector<int> &dimensions,
                                       const std::vector<std::vector<int>> &scattering) const {
    std::vector<int> cells;
    for (auto i = 0u; i < dimensions.size(); ++i) {
      cells.push_back(dimensions[i] == 1 ? VTKWriterUtils<1>::NumberOfCells({scattering[i][0]}) :
                      (dimensions[i] == 2 ? VTKWriterUtils<2>::NumberOfCells({scattering[i][0], scattering[i][1]}) :
                       VTKWriterUtils<3>::NumberOfCells({scattering[i][0], scattering[i][1], scattering[i][2]})));
    }
    return cells;
  }

  int GetNumberOfCellEntries(const std::vector<int> &dimensions, const std::vector<int> &cells) const {
    int sum = 0;
    for (auto i = 0u; i < dimensions.size(); ++i) {
      sum += cells[i] * (dimensions[i] == 1 ? 3 : (dimensions[i] == 2 ? 5 : 9));
    }
    return sum;
  }

  void AddPoints(std::ofstream &file, const std::any &spline, const std::vector<int> &scattering, int spl_dim) const {
    if (spl_dim == 1) {
      VTKWriterUtils<1>::WritePoints(file, spline, {scattering[0]});
    } else if (spl_dim == 2) {
      VTKWriterUtils<2>::WritePoints(file, spline, {scattering[0], scattering[1]});
    } else if (spl_dim == 3) {
      VTKWriterUtils<3>::WritePoints(file, spline, {scattering[0], scattering[1], scattering[2]});
    } else {
      throw std::runtime_error("Only splines of dimensions 1 to 3 can be written to a vtk file.");
    }
  }

  void AddCells(std::ofstream &file, const std::vector<int> &scattering, int offset, int spl_dim) const {
    if (spl_dim == 1) {
      Write1DCells(file, scattering[0], offset);
    } else if (spl_dim == 2) {
      Write2DCells(file, {scattering[0], scattering[1]}, offset);
    } else if (spl_dim == 3) {
      Write3DCells(file, {scattering[0], scattering[1], scattering[2]}, offset);
    }
  }

  void AddCellTypes(std::ofstream &file, const std::vector<int> &scattering, int spl_dim) const {
    if (spl_dim == 1) {
      VTKWriterUtils<1>::WriteCellTypes(file, {scattering[0]}, 3);
    } else if (spl_dim == 2) {
      VTKWriterUtils<2>::WriteCellTypes(file, {scattering[0], scattering[1]}, 9);
    } else if (spl_dim == 3) {
      VTKWriterUtils<3>::WriteCellTypes(file, {scattering[0], scattering[1], scattering[2]}, 12);
    }
  }

  void Write1DCells(std::ofstream &file, int scattering, int offset) const {
    for (int i = 0; i < scattering; ++i) {
      file << "2 " << i + offset << " " << i + 1 + offset << "\n";
    }
  }

  void Write2DCells(std::ofstream &file, std::array<int, 2> scattering, int offset) const {
    util::MultiIndexHandler<2> point_handler({scattering[0] + 1, scattering[1] + 1});
    for (int j = 0; j < scattering[0]; ++j) {
      for (int i = 0; i < scattering[1]; ++i) {
        file << "4 " << (scattering[0] + 1) * j + i + offset << " " << (scattering[0] + 1) * j + i + 1 + offset << " "
             << (scattering[0] + 1) * (j + 1) + i + 1 + offset << " " << (scattering[0] + 1) * (j + 1) + i + offset
             << "\n";
      }
    }
  }

  void Write3DCells(std::ofstream &file, std::array<int, 3> scattering, int offset) const {
    util::MultiIndexHandler<3> point_handler({scattering[0] + 1, scattering[1] + 1, scattering[2] + 1});
    for (; point_handler.Get1DIndex() < point_handler.Get1DLength() - 1; ++point_handler) {
      if (point_handler.GetIndices()[0] != scattering[0] && point_handler.GetIndices()[1] != scattering[1]
          && point_handler.GetIndices()[2] != scattering[2]) {
        file << "8 " << point_handler.Get1DIndex() + offset << " " << (point_handler + 1).Get1DIndex() + offset << " "
             << (point_handler + scattering[0] + 1).Get1DIndex() + offset << " "
             << (point_handler - 1).Get1DIndex() + offset << " "
             << (point_handler + scattering[1] * (scattering[0] + 1)).Get1DIndex() + offset << " "
             << (point_handler + 1).Get1DIndex() + offset << " "
             << (point_handler + scattering[0] + 1).Get1DIndex() + offset << " "
             << (point_handler - 1).Get1DIndex() + offset << std::endl;
        point_handler - (scattering[0] + 1) * (scattering[1] + 2);
      }
    }
  }
};
}  // namespace io

#endif  // SRC_IO_VTK_WRITER_H_
