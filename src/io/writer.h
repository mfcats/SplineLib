/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_WRITER_H_
#define SRC_IO_WRITER_H_

#include <any>
#include <vector>

namespace io {
class Writer {
 public:
  Writer() = default;

  void WriteFile(const std::vector<std::any> &splines,
                 const char *filename,
                 const std::vector<std::vector<int>> &scattering = {}) const;

  std::vector<std::any> GetSplinesOfCorrectDimension(const std::vector<std::any> &splines, int max_dim) const;
  std::vector<int> GetSplinePositionsOfCorrectDimension(const std::vector<std::any> &splines, int max_dim) const;

  void PrintWarningForOmittedSplines(size_t splines, size_t count, int max_dim, const char *filename) const;
};
}  // namespace io

#endif  // SRC_IO_WRITER_H_
