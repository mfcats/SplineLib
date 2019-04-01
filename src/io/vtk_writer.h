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
#include <array>
#include <fstream>
#include <vector>

namespace io {
class VTKWriter {
 public:
  VTKWriter() = default;

  void WriteFile(const std::vector<std::any> &splines,
                 const std::string &filename,
                 const std::vector<std::vector<int>> &scattering) const;

 private:
  std::vector<int> GetSplineDimensions(const std::vector<std::any> &splines) const;

  void ThrowIfScatteringHasWrongSizes(std::vector<std::vector<int>> scattering, std::vector<int> dimensions) const;

  std::vector<int> GetNumberOfAllPoints(const std::vector<int> &dimensions,
                                        const std::vector<std::vector<int>> &scattering) const;
  std::vector<int> GetNumberOfAllCells(const std::vector<int> &dimensions,
                                       const std::vector<std::vector<int>> &scattering) const;
  int GetNumberOfCellEntries(const std::vector<int> &dimensions, const std::vector<int> &cells) const;

  void AddPoints(std::ofstream &file, const std::any &spline, const std::vector<int> &scattering, int spl_dim) const;
  void AddCells(std::ofstream &file, const std::vector<int> &scattering, int offset, int spl_dim) const;
  void AddCellTypes(std::ofstream &file, const std::vector<int> &scattering, int spl_dim) const;

  void Write1DCells(std::ofstream &file, int scattering, int offset) const;
  void Write2DCells(std::ofstream &file, std::array<int, 2> scattering, int offset) const;
  void Write3DCells(std::ofstream &file, std::array<int, 3> scattering, int offset) const;
};
}  // namespace io

#endif  // SRC_IO_VTK_WRITER_H_
