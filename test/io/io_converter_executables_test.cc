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

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

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

std::string GetPathToInstallDir() {
  char exePath[PATH_MAX];
  uint32_t len = sizeof(exePath);
#ifdef __linux__
  char szTmp[32];
  snprintf(szTmp, 32, "/proc/%d/exe", getpid());
  int bytes = readlink(szTmp, exePath, len);
  if (bytes >= 0) {
    exePath[bytes] = '\0';
  }
#elif __APPLE__
  if (_NSGetExecutablePath(exePath, &len) != 0) {
    exePath[0] = '\0';
  } else {
    char *canonicalPath = realpath(exePath, nullptr);
    if (canonicalPath != nullptr) {
      strncpy(exePath, canonicalPath, len);
      free(canonicalPath);
    }
  }
#else
  throw std::runtime_error("Cannot find path of executable.");
#endif
  std::string path(exePath);
  return path.substr(0, path.length() - 14) + "/../src/io/";
}

void CreateLogFile(const char *executable, const char *output, const char *options = "all") {
  std::ofstream outfile("log.txt");
  outfile << "input:\n" << executable << "\n\noutput:\n" << output << "\n\noptions:\n" << options;
  outfile.close();
}

std::string GetLog() {
  std::ifstream newFile("log.txt");
  std::string line, file;
  while (getline(newFile, line)) {
    file += line + "\n";
  }
  return file;
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

 protected:
  io::IGESReader iges_reader_;
  io::IRITReader irit_reader_;
};

TEST_F(Iges2iritExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "iges2irit -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "iges2irit --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

TEST_F(Iges2iritExecutable, Works) {  // NOLINT
  CreateLogFile(iges_read, "out.itd");
  std::vector<std::any> iges_splines = iges_reader_.ReadFile(iges_read);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[0]);
  auto iges_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines[1]);
  std::system((GetPathToInstallDir() + "iges2irit log.txt").c_str());
  std::vector<std::any> irit_splines = irit_reader_.ReadFile("out.itd");
  auto irit_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[0]);
  auto irit_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(irit_splines[1]);
  ASSERT_THAT(irit_nurbs_2d->AreEqual(*iges_nurbs_2d.get(), 1e-6), true);
  ASSERT_THAT(irit_bspline_1d->AreEqual(*iges_bspline_1d.get(), 1e-6), true);
  std::string log = GetLog();
  ASSERT_THAT(log.find("log:\n"), Ne(std::string::npos));
  ASSERT_THAT(
      log.find(std::string("The splines at positions 0, 1 in file ") + iges_read + " have been written to out.itd"),
      Ne(std::string::npos));
  remove("log.txt");
  remove("out.itd");
}

class Iges2vtkExecutable : public Test {
 public:
  Iges2vtkExecutable() = default;
};

TEST_F(Iges2vtkExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "iges2vtk -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "iges2vtk --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

class Iges2xmlExecutable : public Test {
 public:
  Iges2xmlExecutable() = default;
};

TEST_F(Iges2xmlExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "iges2xml -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "iges2xml --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

class Irit2igesExecutable : public Test {
 public:
  Irit2igesExecutable() = default;
};

TEST_F(Irit2igesExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "irit2iges -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "irit2iges --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

class Irit2vtkExecutable : public Test {
 public:
  Irit2vtkExecutable() = default;
};

TEST_F(Irit2vtkExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "irit2vtk -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "irit2vtk --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

class Irit2xmlExecutable : public Test {
 public:
  Irit2xmlExecutable() = default;
};

TEST_F(Irit2xmlExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "irit2xml -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "irit2xml --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

class Xml2igesExecutable : public Test {
 public:
  Xml2igesExecutable() = default;
};

TEST_F(Xml2igesExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "xml2iges -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "xml2iges --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

class Xml2iritExecutable : public Test {
 public:
  Xml2iritExecutable() = default;
};

TEST_F(Xml2iritExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "xml2irit -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "xml2irit --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

class Xml2vtkExecutable : public Test {
 public:
  Xml2vtkExecutable() = default;
};

TEST_F(Xml2vtkExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "xml2vtk -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "xml2vtk --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}
