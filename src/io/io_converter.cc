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

void io::IOConverter::ConvertFile(const char *input_filename,
                                  const char *output_filename,
                                  const std::vector<int> &positions,
                                  const std::vector<std::vector<int>> &scattering) {
  std::vector<std::any> splines = ReadFile(input_filename);
  WriteFile(splines, output_filename, scattering);
}

std::vector<std::any> io::IOConverter::ReadFile(const char *filename) {
  if (util::StringOperations::EndsWith(filename, ".iges")) {
    io::IGESReader iges_reader;
    return iges_reader.ReadFile(filename);
  } else if (util::StringOperations::EndsWith(filename, ".itd")) {
    io::IRITReader irit_reader;
    return irit_reader.ReadFile(filename);
  } else if (util::StringOperations::EndsWith(filename, ".xml")) {
    io::XMLReader xml_reader;
    return xml_reader.ReadFile(filename);
  } else {
    throw std::runtime_error(R"(Only files of format ".iges", ".itd" and ".xml" can be read.)");
  }
}

void io::IOConverter::WriteFile(const std::vector<std::any> &splines,
                                const char *filename,
                                const std::vector<std::vector<int>> &scattering) const {
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

std::vector<std::any> io::IOConverter::GetSplinesOfCorrectDimension(const std::vector<std::any> &splines,
                                                                    int max_dim) const {
  std::vector<int> positions = GetSplinePositionsOfCorrectDimension(splines, max_dim);
  std::vector<std::any> splines_with_max_dim = util::VectorUtils<std::any>::FilterVector(splines, positions);
  return splines_with_max_dim;
}

std::vector<int> io::IOConverter::GetSplinePositionsOfCorrectDimension(const std::vector<std::any> &splines,
                                                                       int max_dim) const {
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
