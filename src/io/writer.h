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
#include "xml_writer.h"

namespace io {
class Writer {
 public:
  Writer() = default;

  virtual void WriteFile(std::vector<std::any> &splines, const char *filename) {
    if (util::StringOperations::EndsWith(filename, ".iges")) {
      io::IGESWriter iges_writer;
      splines = GetSplinesOfCorrectDimensions(splines, 2);
      iges_writer.WriteFile(splines, filename);
    } else if (util::StringOperations::EndsWith(filename, ".itd")) {
      io::IRITWriter irit_writer;
      splines = GetSplinesOfCorrectDimensions(splines, 3);
      return irit_writer.WriteFile(splines, filename);
    } else if (util::StringOperations::EndsWith(filename, ".xml")) {
      io::XMLWriter xml_writer;
      splines = GetSplinesOfCorrectDimensions(splines, 4);
      return xml_writer.WriteFile(splines, filename);
    } else {
      throw std::runtime_error("Only files of format iges, itd and xml can be written.");
    }
  }

 private:
  std::vector<std::any> GetSplinesOfCorrectDimensions(const std::vector<std::any> &splines, int max_dim) {
    std::vector<std::any> splines_with_max_dim;
    for (const auto &spline : splines) {
      if (util::AnyCasts::GetSplineDimension(spline) <= max_dim) {
        splines_with_max_dim.emplace_back(spline);
      }
    }
    return splines_with_max_dim;
  }
};
}  // namespace io

#endif  // SRC_IO_WRITER_H_
