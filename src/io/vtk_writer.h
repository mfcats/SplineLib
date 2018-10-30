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
      newFile << "# vtk DataFile Version 3.0" << std::endl << "Spline from Splinelib\n" << "ASCII\n\n";
      for (auto i = 0u; i < splines.size(); ++i) {
        AddSpline(newFile, splines[i], scattering[i]);
      }
      newFile.close();
    }
  }

 private:
  void AddSpline(std::ofstream &file, const std::any &spline, const std::vector<int> &scattering) const {
    int spline_dimension = util::AnyCasts::GetSplineDimension(spline);
    switch (spline_dimension) {
      case 1: {
        Write1DSpline(file, spline, scattering[0]);
        break;
      }
      case 2: {
        Write2DSpline(file, spline, {scattering[0], scattering[1]});
        break;
      }
      case 3: {
        Write3DSpline(file, spline, {scattering[0], scattering[1], scattering[2]});
        break;
      }
      default: {
        throw std::runtime_error("Only splines of dimensions 1 to 3 can be written to a vtk file.");
      }
    }
  }

  void Write1DSpline(std::ofstream &file, const std::any &spline, int scattering) const {
    std::shared_ptr<spl::Spline<1>> spline_ptr = util::AnyCasts::GetSpline<1>(spline);
    double lowest_knot = spline_ptr->GetKnotVector(0)->GetKnot(0).get();
    double highest_knot = spline_ptr->GetKnotVector(0)->GetLastKnot().get();
    file << "DATASET POLYDATA\nPOINTS " << scattering + 1 << " double\n";
    for (int i = 0; i <= scattering; ++i) {
      std::array<ParamCoord, 1> param_coord = {ParamCoord(lowest_knot + i * (highest_knot - lowest_knot) / scattering)};
      for (int j = 0; j < 3; ++j) {
        if (j < spline_ptr->GetDimension()) file << spline_ptr->Evaluate(param_coord, {j})[0] << " ";
        else file << 0 << " ";
      }
      file << "\n";
    }
    file << "\nLINES " << scattering << " " << 3 * scattering << "\n";
    for (int i = 0; i < scattering; ++i) {
      file << "2 " << i << " " << i + 1 << "\n";
    }
  }

  void Write2DSpline(std::ofstream &file, const std::any &spline, std::array<int, 2> scattering) const {
    std::shared_ptr<spl::Spline<2>> spline_ptr = util::AnyCasts::GetSpline<2>(spline);
    std::array<double, 2> lowest_knot = {spline_ptr->GetKnotVector(0)->GetKnot(0).get(),
                                         spline_ptr->GetKnotVector(1)->GetKnot(0).get()};
    std::array<double, 2> highest_knot = {spline_ptr->GetKnotVector(0)->GetLastKnot().get(),
                                          spline_ptr->GetKnotVector(1)->GetLastKnot().get()};
    file << "DATASET POLYDATA\nPOINTS " << (scattering[0] + 1) * (scattering[1] + 1) << " double\n";
    for (int i = 0; i <= scattering[1]; ++i) {
      for (int j = 0; j <= scattering[0]; ++j) {
        std::array<ParamCoord, 2> param_coord =
            {ParamCoord(lowest_knot[0] + j * (highest_knot[0] - lowest_knot[0]) / scattering[0]),
             ParamCoord(lowest_knot[1] + i * (highest_knot[1] - lowest_knot[1]) / scattering[1])};
        for (int k = 0; k < 3; ++k) {
          if (k < spline_ptr->GetDimension()) file << spline_ptr->Evaluate(param_coord, {k})[0] << " ";
          else file << 0 << " ";
        }
        file << "\n";
      }
    }
    file << "\nPOLYGONS " << scattering[0] * scattering[1] << " " << 5 * scattering[0] * scattering[1] << "\n";
    for (int j = 0; j < scattering[0]; ++j) {
      for (int i = 0; i < scattering[1]; ++i) {
        file << "4 " << (scattering[0] + 1) * j + i << " " << (scattering[0] + 1) * j + i + 1 << " " <<
             (scattering[0] + 1) * (j + 1) + i + 1 << " " << (scattering[0] + 1) * (j + 1) + i << "\n";
      }
    }
  }

  void Write3DSpline(std::ofstream &file, const std::any &spline, std::array<int, 3> scattering) const {
    std::shared_ptr<spl::Spline<3>> spline_ptr = util::AnyCasts::GetSpline<3>(spline);
    std::array<double, 3> lowest_knot = {spline_ptr->GetKnotVector(0)->GetKnot(0).get(),
                                         spline_ptr->GetKnotVector(1)->GetKnot(0).get(),
                                         spline_ptr->GetKnotVector(2)->GetKnot(0).get()};
    std::array<double, 3> highest_knot = {spline_ptr->GetKnotVector(0)->GetLastKnot().get(),
                                          spline_ptr->GetKnotVector(1)->GetLastKnot().get(),
                                          spline_ptr->GetKnotVector(2)->GetLastKnot().get()};
    file << "DATASET UNSTRUCTURED_GRID\nPOINTS " << (scattering[0] + 1) * (scattering[1] + 1) * (scattering[2] + 1)
         << " double\n";
    for (int i = 0; i <= scattering[2]; ++i) {
      for (int j = 0; j <= scattering[1]; ++j) {
        for (int k = 0; k <= scattering[0]; ++k) {
          std::array<ParamCoord, 3> param_coord =
              {ParamCoord(lowest_knot[0] + k * (highest_knot[0] - lowest_knot[0]) / scattering[0]),
               ParamCoord(lowest_knot[1] + j * (highest_knot[1] - lowest_knot[1]) / scattering[1]),
               ParamCoord(lowest_knot[2] + i * (highest_knot[2] - lowest_knot[2]) / scattering[2])};
          for (int l = 0; l < 3; ++l) {
            if (l < spline_ptr->GetDimension()) file << spline_ptr->Evaluate(param_coord, {l})[0] << " ";
            else file << 0 << " ";
          }
          file << "\n";
        }
      }
    }
    file << "\nCELLS " << scattering[0] * scattering[1] * scattering[2] << " "
         << 9 * scattering[0] * scattering[1] * scattering[2] << "\n";
    util::MultiIndexHandler<3> point_handler({scattering[0] + 1, scattering[1] + 1, scattering[2] + 1});
    for (point_handler.SetIndices({0, 0, 0}); point_handler.Get1DIndex() < point_handler.Get1DLength() - 1;
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
    file << "\nCELL_TYPES " << scattering[0] * scattering[1] * scattering[2] << "\n";
    for (int j = 0; j < scattering[0] * scattering[1] * scattering[2]; ++j) {
      file << "12\n";
    }
  }

};
}  // namespace io
#endif  // SRC_IO_VTK_WRITER_H_
