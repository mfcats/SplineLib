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

#include <vector>

#include "iges_writer.h"
#include "irit_writer.h"
#include "string_operations.h"
#include "vtk_writer.h"
#include "xml_writer.h"

namespace io {
class Writer {
 public:
  Writer() = default;

  void WriteFile(const std::vector<std::any> &splines,
                 const char *filename,
                 const std::vector<std::vector<int>> &scattering = {}) const {
    if (util::StringOperations::EndsWith(filename, ".iges")) {
      io::IGESWriter iges_writer;
      std::vector<std::any> splines_with_max_dim = GetSplinesOfCorrectDimension(splines, 2);
      iges_writer.WriteFile(splines_with_max_dim, filename);
      PrintWarningForOmittedSplines(splines.size(), splines_with_max_dim.size(), 2, filename);
    } else if (util::StringOperations::EndsWith(filename, ".itd")) {
      io::IRITWriter irit_writer;
      std::vector<std::any> splines_with_max_dim = GetSplinesOfCorrectDimension(splines, 3);
      PrintWarningForOmittedSplines(splines.size(), splines_with_max_dim.size(), 3, filename);
      return irit_writer.WriteFile(splines_with_max_dim, filename);
    } else if (util::StringOperations::EndsWith(filename, ".vtk")) {
      io::VTKWriter vtk_writer;
      std::vector<std::any> splines_with_max_dim = GetSplinesOfCorrectDimension(splines, 3);
      PrintWarningForOmittedSplines(splines.size(), splines_with_max_dim.size(), 3, filename);
      return vtk_writer.WriteFile(splines_with_max_dim, filename, scattering);
    } else if (util::StringOperations::EndsWith(filename, ".xml")) {
      io::XMLWriter xml_writer;
      std::vector<std::any> splines_with_max_dim = GetSplinesOfCorrectDimension(splines, 4);
      PrintWarningForOmittedSplines(splines.size(), splines_with_max_dim.size(), 4, filename);
      return xml_writer.WriteFile(splines_with_max_dim, filename);
    } else {
      throw std::runtime_error(R"(Only files of format ".iges", ".itd", ".vtk" and ".xml" can be written.)");
    }
  }

 private:
  std::vector<std::any> GetSplinesOfCorrectDimension(const std::vector<std::any> &splines, int max_dim) const {
    std::vector<std::any> splines_with_max_dim;
    for (const auto &spline : splines) {
      if (util::AnyCasts::GetSplineDimension(spline) <= max_dim) {
        splines_with_max_dim.emplace_back(spline);
      }
    }
    return splines_with_max_dim;
  }

  void PrintWarningForOmittedSplines(size_t splines, size_t count, int max_dim, const char *filename) const {
    if (count < splines) {
      std::cerr << "Only the " << count << " splines of dimension equal or less than " << max_dim
                << " have been written to " << filename << "." << std::endl;
    }
  }
};
}  // namespace io

#endif  // SRC_IO_WRITER_H_
