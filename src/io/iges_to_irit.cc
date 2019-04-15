/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "converter_log.h"
#include "iges_reader.h"
#include "writer.h"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    throw std::runtime_error("Exactly one name of the log file is required as command line argument.");
  }
  const char *log_file = argv[1];  // NOLINT
  if (util::StringOperations::StartsWith(log_file, "--help") || util::StringOperations::StartsWith(log_file, "-h")) {
    io::ConverterLog help;
    help.PrintHelp();
    return 0;
  }
  io::ConverterLog log(log_file);
  std::vector<std::any> splines;

  try {
    io::IGESReader iges_reader;
    splines = iges_reader.ReadFile(log.GetInput());
  } catch (std::runtime_error &error) {
    std::cerr << error.what() << std::endl;
    throw error;
  } catch (...) {
    throw std::runtime_error(R"(The input file isn't of correct ".iges" format.)");
  }

  io::IRITWriter irit_writer;
  io::Writer writer;
  std::vector<int> positions = log.GetPositions(writer.GetSplinePositionsOfCorrectDimension(splines, 3));
  std::vector<std::any> splines_with_max_dim = util::VectorUtils<std::any>::FilterVector(splines, positions);
  irit_writer.WriteFile(splines_with_max_dim, log.GetOutput());

  log.WriteLog();
  return 0;
}