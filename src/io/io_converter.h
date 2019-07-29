/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_IO_CONVERTER_H_
#define SRC_IO_IO_CONVERTER_H_

#include <any>
#include <string>
#include <vector>

#include "any_casts.h"
#include "iges_reader.h"
#include "iges_writer.h"
#include "irit_reader.h"
#include "irit_writer.h"
#include "string_operations.h"
#include "vtk_writer.h"
#include "xml_reader.h"
#include "xml_writer.h"

namespace io {
class IOConverter {
 public:
  IOConverter(const char *input_filename, const char *output_filename);

  std::vector<int> ConvertFile(const std::vector<int> &positions = {},
                               const std::vector<std::vector<int>> &scattering = {});

  static std::vector<int> GetSplinePositionsOfCorrectDimension(const std::vector<std::any> &splines, int max_dim);

  enum file_format { error, iges, irit, vtk, xml };

 private:
  file_format GetFileFormat(const char *filename) const;

  io::Reader *GetReader() const;
  void GetWriter(const std::vector<std::any> &splines, const std::vector<int> &positions,
                 const std::vector<std::vector<int>> &scattering);
  void WriteFile(const std::vector<std::any> &splines, const std::vector<int> &positions, int max_dim,
                 const std::shared_ptr<io::Writer> &writer);

  void GetPositions(const std::vector<int> &positions, const std::vector<int> &possible_positions);
  std::vector<std::vector<int>> GetScattering(const std::vector<std::vector<int>> &scattering,
                                              const std::vector<int> &positions,
                                              const std::vector<int> &written);

  const char *input_filename_;
  const char *output_filename_;
  file_format input_format_;
  file_format output_format_;
  std::vector<int> written_;
};
}  // namespace io

#endif  // SRC_IO_IO_CONVERTER_H_
