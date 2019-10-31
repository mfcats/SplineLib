/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "io_converter.h"

namespace splinelib::src::io {
IOConverter::IOConverter(const char *input_filename, const char *output_filename)
    : input_filename_(input_filename),
      output_filename_(output_filename),
      input_format_(GetFileFormat(input_filename)),
      output_format_(GetFileFormat(output_filename)) {}

std::vector<int> IOConverter::ConvertFile(const std::vector<int> &positions,
                                              const std::vector<std::vector<int>> &scattering) {
  std::vector<std::any> splines = GetReader()->ReadFile(input_filename_);
  GetWriter(splines, positions, scattering);
  return written_;
}

std::vector<int> IOConverter::GetSplinePositionsOfCorrectDimension(const std::vector<std::any> &splines,
                                                                       int max_dim) {
  std::vector<int> spline_positions_with_max_dim;
  for (size_t i = 0; i < splines.size(); ++i) {
    if (util::AnyCasts::GetSplineDimension(splines[i]) <= max_dim) {
      spline_positions_with_max_dim.emplace_back(i);
    }
  }
  return spline_positions_with_max_dim;
}

IOConverter::file_format IOConverter::GetFileFormat(const char *filename) const {
  if (util::string_operations::EndsWith(filename, ".iges")) {
    return iges;
  }
  if (util::string_operations::EndsWith(filename, ".itd")) {
    return irit;
  }
  if (util::string_operations::EndsWith(filename, ".vtk")) {
    return vtk;
  }
  if (util::string_operations::EndsWith(filename, ".xml")) {
    return xml;
  }
  return error;
}

Reader *IOConverter::GetReader() const {
  if (input_format_ == iges) {
    return new IGESReader;
  }
  if (input_format_ == irit) {
    return new IRITReader;
  }
  if (input_format_ == xml) {
    return new XMLReader;
  }
  throw std::runtime_error(R"(Only files of format ".iges", ".itd" and ".xml" can be read.)");
}

void IOConverter::GetWriter(const std::vector<std::any> &splines, const std::vector<int> &positions,
                                const std::vector<std::vector<int>> &scattering) {
  if (output_format_ == iges) {
    std::shared_ptr<IGESWriter> ptr = std::make_shared<IGESWriter>(IGESWriter());
    WriteFile(splines, positions, 2, std::dynamic_pointer_cast<Writer>(ptr));
  } else if (output_format_ == irit) {
    std::shared_ptr<IRITWriter> ptr = std::make_shared<IRITWriter>(IRITWriter());
    WriteFile(splines, positions, 3, std::dynamic_pointer_cast<Writer>(ptr));
  } else if (output_format_ == vtk) {
    std::shared_ptr<VTKWriter> ptr = std::make_shared<VTKWriter>(VTKWriter());
    GetPositions(positions, GetSplinePositionsOfCorrectDimension(splines, 3));
    std::vector<std::any> splines_with_max_dim = util::vector_utils::GetEntriesAtIndices<std::any>(splines, written_);
    ptr->WriteFile(splines_with_max_dim, output_filename_, GetScattering(scattering, positions, written_));
  } else if (output_format_ == xml) {
    std::shared_ptr<XMLWriter> ptr = std::make_shared<XMLWriter>(XMLWriter());
    WriteFile(splines, positions, 4, std::dynamic_pointer_cast<Writer>(ptr));
  } else {
    throw std::runtime_error(R"(Only files of format ".iges", ".itd", ".vtk" and ".xml" can be written.)");
  }
}

void IOConverter::WriteFile(const std::vector<std::any> &splines, const std::vector<int> &positions, int max_dim,
                                const std::shared_ptr<Writer> &writer) {
  GetPositions(positions, GetSplinePositionsOfCorrectDimension(splines, max_dim));
  std::vector<std::any> splines_with_max_dim = util::vector_utils::GetEntriesAtIndices<std::any>(splines, written_);
  writer->WriteFile(splines_with_max_dim, output_filename_);
}

void IOConverter::GetPositions(const std::vector<int> &positions,
                                   const std::vector<int> &possible_positions) {
  if (positions.empty()) {
    written_ = possible_positions;
  } else {
    for (const auto &pos : positions) {
      if (std::find(possible_positions.begin(), possible_positions.end(), pos) != possible_positions.end()) {
        written_.emplace_back(pos);
      }
    }
  }
}

std::vector<std::vector<int>> IOConverter::GetScattering(const std::vector<std::vector<int>> &scattering,
                                                         const std::vector<int> &positions,
                                                         const std::vector<int> &written) {
  std::vector<int> indices;
  if (!positions.empty()) {
    for (int i = 0; i < static_cast<int>(positions.size()); ++i) {
      if (std::find(written.begin(), written.end(), positions[i]) != written.end()) {
        indices.emplace_back(i);
      }
    }
  } else {
    indices = written;
  }
  return util::vector_utils::GetEntriesAtIndices(scattering, indices);
}
}  // namespace splinelib::src::io
