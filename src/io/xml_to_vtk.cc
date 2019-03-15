/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "xml_reader.h"
#include "writer.h"

int main(int argc, char *argv[]) {
  if (argc != 3) {
    throw std::runtime_error("Exactly one name of the input file and of the output file are required");
  }
  const char *input = argv[1];
  const char *output = argv[2];

  std::vector<std::any> splines;
  try {
    io::XMLReader xml_reader;
    splines = xml_reader.ReadFile(input);
  } catch (...) {
    throw std::runtime_error(R"(The input file isn't of correct ".xml" format.)");
  }
  io::VTKWriter vtk_writer;
  io::Writer writer;
  std::vector<std::any> splines_with_max_dim = writer.GetSplinesOfCorrectDimension(splines, 3);
  vtk_writer.WriteFile(splines_with_max_dim, output, {});
  writer.PrintWarningForOmittedSplines(splines.size(), splines_with_max_dim.size(), 3, output);
  return 0;
}
