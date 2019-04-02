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

io::IOConverter::IOConverter(const char *input_filename, const char *output_filename)
    : input_filename_(input_filename),
      output_filename_(output_filename),
      input_format_(GetFileFormat(input_filename)),
      output_format_(GetFileFormat(output_filename)) {}

void io::IOConverter::ConvertFile(const std::vector<int> &positions,
                                  const std::vector<std::vector<int>> &scattering) const {
  std::vector<std::any> splines = ReadFile();
  WriteFile(splines, scattering);
}

std::vector<std::any> io::IOConverter::ReadFile() const {
  if (input_format_ == iges) {
    io::IGESReader iges_reader;
    return iges_reader.ReadFile(input_filename_);
  } else if (input_format_ == irit) {
    io::IRITReader irit_reader;
    return irit_reader.ReadFile(input_filename_);
  } else if (input_format_ == xml) {
    io::XMLReader xml_reader;
    return xml_reader.ReadFile(input_filename_);
  } else {
    throw std::runtime_error(R"(Only files of format ".iges", ".itd" and ".xml" can be read.)");
  }
}

void io::IOConverter::WriteFile(const std::vector<std::any> &splines,
                                const std::vector<std::vector<int>> &scattering) const {
  if (output_format_ == iges) {
    io::IGESWriter iges_writer;
    std::vector<std::any> splines_with_max_dim = GetSplinesOfCorrectDimension(splines, 2);
    iges_writer.WriteFile(splines_with_max_dim, output_filename_);
    PrintWarningForOmittedSplines(splines.size(), splines_with_max_dim.size(), 2, output_filename_);
  } else if (output_format_ == irit) {
    io::IRITWriter irit_writer;
    std::vector<std::any> splines_with_max_dim = GetSplinesOfCorrectDimension(splines, 3);
    PrintWarningForOmittedSplines(splines.size(), splines_with_max_dim.size(), 3, output_filename_);
    return irit_writer.WriteFile(splines_with_max_dim, output_filename_);
  } else if (output_format_ == vtk) {
    io::VTKWriter vtk_writer;
    std::vector<std::any> splines_with_max_dim = GetSplinesOfCorrectDimension(splines, 3);
    PrintWarningForOmittedSplines(splines.size(), splines_with_max_dim.size(), 3, output_filename_);
    return vtk_writer.WriteFile(splines_with_max_dim, output_filename_, scattering);
  } else if (output_format_ == xml) {
    io::XMLWriter xml_writer;
    std::vector<std::any> splines_with_max_dim = GetSplinesOfCorrectDimension(splines, 4);
    PrintWarningForOmittedSplines(splines.size(), splines_with_max_dim.size(), 4, output_filename_);
    return xml_writer.WriteFile(splines_with_max_dim, output_filename_);
  } else {
    throw std::runtime_error(R"(Only files of format ".iges", ".itd", ".vtk" and ".xml" can be written.)");
  }
}

std::vector<std::any> io::IOConverter::GetSplinesOfCorrectDimension(const std::vector<std::any> &splines,
                                                                    int max_dim) const {
  std::vector<int> positions = GetSplinePositionsOfCorrectDimension(splines, max_dim);
  std::vector<std::any> splines_with_max_dim = util::VectorUtils<std::any>::FilterVector(splines, positions);
  return splines_with_max_dim;
}

std::vector<int> io::IOConverter::GetSplinePositionsOfCorrectDimension(const std::vector<std::any> &splines,
                                                                       int max_dim) {
  std::vector<int> spline_positions_with_max_dim;
  for (size_t i = 0; i < splines.size(); ++i) {
    if (util::AnyCasts::GetSplineDimension(splines[i]) <= max_dim) {
      spline_positions_with_max_dim.emplace_back(i);
    }
  }
  return spline_positions_with_max_dim;
}

void io::IOConverter::PrintWarningForOmittedSplines(size_t splines,
                                                    size_t count,
                                                    int max_dim,
                                                    const char *filename) const {
  if (count < splines) {
    std::cerr << "Only the " << count << " splines of dimension equal or less than " << max_dim
              << " have been written to " << filename << "." << std::endl;
  }
}

std::vector<int> io::IOConverter::GetPositions(const std::vector<int> &positions,
                                               const std::vector<int> &possible_positions) {
  std::vector<int> written;
  if (positions.size() == 0) {
    written = possible_positions;
  } else {
    for (const auto &pos : positions) {
      if (std::find(possible_positions.begin(), possible_positions.end(), pos) != possible_positions.end()) {
        written.emplace_back(pos);
      } else {
        not_written_.emplace_back(pos);
      }
    }
  }
  return written;
}

io::IOConverter::file_format io::IOConverter::GetFileFormat(const char *filename) const {
  if (util::StringOperations::EndsWith(filename, ".iges")) {
    return iges;
  } else if (util::StringOperations::EndsWith(filename, ".itd")) {
    return irit;
  } else if (util::StringOperations::EndsWith(filename, ".vtk")) {
    return vtk;
  } else if (util::StringOperations::EndsWith(filename, ".xml")) {
    return xml;
  } else {
    return error;
  }
}
