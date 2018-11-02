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
      std::vector<int> dimensions = GetSplineDimensions(splines);
      std::vector<int> cells = GetNumberOfAllCells(dimensions, scattering);
      std::vector<int> points = GetNumberOfAllPoints(dimensions, scattering);
      newFile << "POINTS " << std::accumulate(points.begin(), points.end(), 0) << " double\n";
      for (auto i = 0u; i < splines.size(); ++i) {
        AddPoints(newFile, splines[i], scattering[i]);
      }
      newFile << "\nCELLS " << std::accumulate(cells.begin(), cells.end(), 0) << " "
              << GetCellEntries(dimensions, cells) << "\n";
      for (auto i = 0u; i < splines.size(); ++i) {
        AddCells(newFile, splines[i], scattering[i], std::accumulate(points.begin(), points.begin() + i, 0));
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

  std::vector<int> GetNumberOfAllPoints(std::vector<int> &dimensions,
                                        const std::vector<std::vector<int>> &scattering) const {
    std::vector<int> points;
    for (auto i = 0u; i < dimensions.size(); ++i) {
      points.push_back(dimensions[i] == 1 ? NumberOfCells<1>({scattering[i][0] + 1}) :
                       (dimensions[i] == 2 ? NumberOfCells<2>({scattering[i][0] + 1, scattering[i][1] + 1}) :
                        NumberOfCells<3>({scattering[i][0] + 1, scattering[i][1] + 1, scattering[i][2] + 1})));
    }
    return points;
  }

  std::vector<int> GetNumberOfAllCells(std::vector<int> &dimensions,
                                       const std::vector<std::vector<int>> &scattering) const {
    std::vector<int> cells;
    for (auto i = 0u; i < dimensions.size(); ++i) {
      cells.push_back(dimensions[i] == 1 ? NumberOfCells<1>({scattering[i][0]}) :
                      (dimensions[i] == 2 ? NumberOfCells<2>({scattering[i][0], scattering[i][1]}) :
                       NumberOfCells<3>({scattering[i][0], scattering[i][1], scattering[i][2]})));
    }
    return cells;
  }

  int GetCellEntries(std::vector<int> &dimensions, std::vector<int> cells) const {
    int sum = 0;
    for (auto i = 0u; i < dimensions.size(); ++i) {
      sum += cells[i] * (dimensions[i] == 1 ? 3 : (dimensions[i] == 2 ? 5 : 9));
    }
    return sum;
  }

  void AddPoints(std::ofstream &file, const std::any &spline, const std::vector<int> &scattering) const {
    int spline_dimension = util::AnyCasts::GetSplineDimension(spline);
    if (spline_dimension == 1) {
      WritePoints<1>(file, spline, {scattering[0]});
    } else if (spline_dimension == 2) {
      WritePoints<2>(file, spline, {scattering[0], scattering[1]});
    } else if (spline_dimension == 3) {
      WritePoints<3>(file, spline, {scattering[0], scattering[1], scattering[2]});
    } else {
      throw std::runtime_error("Only splines of dimensions 1 to 3 can be written to a vtk file.");
    }
  }

  void AddCells(std::ofstream &file, const std::any &spline, const std::vector<int> &scattering, int offset) const {
    int spline_dimension = util::AnyCasts::GetSplineDimension(spline);
    if (spline_dimension == 1) {
      Write1DCells(file, scattering[0], offset);
    } else if (spline_dimension == 2) {
      Write2DCells(file, {scattering[0], scattering[1]}, offset);
    } else if (spline_dimension == 3) {
      Write3DCells(file, {scattering[0], scattering[1], scattering[2]}, offset);
    }
  }

  void AddCellTypes(std::ofstream &file, const std::vector<int> &scattering, int spline_dimension) const {
    if (spline_dimension == 1) {
      WriteCellTypes<1>(file, {scattering[0]}, 3);
    } else if (spline_dimension == 2) {
      WriteCellTypes<2>(file, {scattering[0], scattering[1]}, 9);
    } else if (spline_dimension == 3) {
      WriteCellTypes<3>(file, {scattering[0], scattering[1], scattering[2]}, 12);
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
    util::MultiIndexHandler<3> point_handler(GetPointHandlerLength<3>(scattering));
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
  int NumberOfCells(std::array<int, DIM> scattering) const {
    int product = 1;
    for (int i = 0; i < DIM; ++i) {
      product *= scattering[i];
    }
    return product;
  }

  template<int DIM>
  std::array<int, DIM> GetPointHandlerLength(std::array<int, DIM> scattering) const {
    for (int i = 0; i < DIM; ++i) {
      ++scattering[i];
    }
    return scattering;
  }

  template<int DIM>
  void WritePoints(std::ofstream &file, const std::any &spline, std::array<int, DIM> scattering) const {
    std::shared_ptr<spl::Spline<DIM>> spline_ptr = util::AnyCasts::GetSpline<DIM>(spline);
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
    for (int i = 0; i < NumberOfCells<DIM>(scattering); ++i) {
      file << cell_type << "\n";
    }
  }
};
}  // namespace io
#endif  // SRC_IO_VTK_WRITER_H_
