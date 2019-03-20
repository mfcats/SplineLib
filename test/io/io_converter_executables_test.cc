/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/


#include <config_iges.h>
#include <config_irit.h>
#include <config_xml.h>
#include <cstdlib>

#include "gmock/gmock.h"

#include "io_converter.h"

using testing::Test;
using testing::Ne;

std::string GetCommandOutput(const std::string &command) {
  std::string result, file;
  FILE *pipe{popen(command.c_str(), "r")};
  char buffer[256];

  while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
    file = buffer;
    result += file.substr(0, file.size() - 1) + "\n";
  }

  pclose(pipe);
  return result;
}

bool CompareToHelpOutput(const std::string &string) {
  std::string help = std::string("The log file has to be of the following format:\n")
      + "input:\n# path to input file\n\n" + "output:\n# path to output file\n\n"
      + "options:\n# list separated by line breaks of positions of splines in input file to be written to the output "
      + "file, 'all' for all splines in input file\n\n"
      + "If this is a converter to VTK format, the log file has to expanded by the following entry:\n" + "scattering:\n"
      + "# scattering for each spline separated by line breaks and for each dimension separated by spaces\n\n"
      + "The log file is expanded by a log entry if converting succeed:\n"
      + "log:\n# time date\n# spline positions in input file that have been written to output file.\n";
  return string == help;
}

class Iges2iritExecutable : public Test {
 public:
  Iges2iritExecutable() = default;
};

TEST_F(Iges2iritExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput("./../src/io/iges2irit -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput("./../src/io/iges2irit --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

TEST_F(Iges2iritExecutable, Works) {  // NOLINT
  std::ofstream outfile("log.txt");
  outfile << "input:\n" << iges_read << "\n\noutput:\nout.itd\n\noptions:\nall";
  outfile.close();
  std::system("./../src/io/iges2irit log.txt");

  std::ifstream newFile("log.txt");
  ASSERT_THAT(newFile.good(), true);
  std::string line, file;
  while (getline(newFile, line)) {
    file += line + "\n";
  }
  ASSERT_THAT(file.find("log:\n"), Ne(std::string::npos));
  ASSERT_THAT(
      file.find(std::string("The splines at positions 0, 1 in file ") + iges_read + " have been written to out.itd"),
      Ne(std::string::npos));
  remove("log.txt");
//  remove("out.itd");
}
