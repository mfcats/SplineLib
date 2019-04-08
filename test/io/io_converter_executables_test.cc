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
  FILE *pipe{popen(command.c_str(), "r")};  // NOLINT
  char buffer[256];
  while (fgets(static_cast<char *>(buffer), sizeof(buffer), pipe) != nullptr) {
    file = static_cast<char *>(buffer);
    result += file.substr(0, file.size() - 1) + "\n";
  }
  pclose(pipe);
  return result;
}

std::string GetPathToInstallDir() {
  char exec_path[PATH_MAX];
  uint32_t length = sizeof(exec_path);
#ifdef __linux__
  char szTmp[32];
  snprintf(static_cast<char*>(szTmp), sizeof(szTmp) / sizeof(char), "/proc/%d/exe", getpid());  // NOLINT
  auto bytes = static_cast<int>(readlink(static_cast<const char*>(szTmp), static_cast<char*>(exec_path), length));
  if (bytes >= 0) {
    exec_path[bytes] = '\0';
  }
#elif __APPLE__
  if (_NSGetExecutablePath(exec_path, &length) != 0) {
    exec_path[0] = '\0';
  } else {
    char *canonical_path = realpath(exec_path, nullptr);
    if (canonical_path != nullptr) {
      strncpy(exec_path, canonical_path, length);
      free(canonical_path);
    }
  }
#else
  throw std::runtime_error("Cannot find path of executable.");
#endif
  std::string path(static_cast<char *>(exec_path));
  return path.substr(0, path.length() - 14) + "/../src/io/";
}

void CreateLogFile(const char *executable, const char *output, const char *options = "all") {
  std::ofstream outfile("log.txt");
  outfile << "input:\n" << executable << "\n\noutput:\n" << output << "\n\noptions:\n" << options;
  outfile.close();
}

std::string GetFileContent(const std::string &filename) {
  std::ifstream newFile(filename);
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

TEST_F(Iges2iritExecutable, WritesLog) {  // NOLINT
  CreateLogFile(iges_read, "out.itd");
  std::system((GetPathToInstallDir() + "iges2irit log.txt").c_str());  // NOLINT
  std::string log = GetFileContent("log.txt");
  ASSERT_THAT(log.find("log:\n"), Ne(std::string::npos));
  ASSERT_THAT(log.find(std::string("The splines at positions 0, 1 in file ") + iges_read
                           + " have been written to out.itd"), Ne(std::string::npos));
  remove("log.txt");
  remove("out.itd");
}

TEST_F(Iges2iritExecutable, ConvertsFileCorrectly) {  // NOLINT
  CreateLogFile(iges_read, "out.itd");
  std::vector<std::any> iges_splines = iges_reader_.ReadFile(iges_read);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[0]);
  auto iges_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines[1]);
  std::system((GetPathToInstallDir() + "iges2irit log.txt").c_str());  // NOLINT
  std::vector<std::any> irit_splines = irit_reader_.ReadFile("out.itd");
  auto irit_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[0]);
  auto irit_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(irit_splines[1]);
  ASSERT_THAT(irit_nurbs_2d->AreEqual(*iges_nurbs_2d.get(), 1e-6), true);
  ASSERT_THAT(irit_bspline_1d->AreEqual(*iges_bspline_1d.get(), 1e-6), true);
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

TEST_F(Iges2vtkExecutable, WritesLog) {  // NOLINT
  CreateLogFile(iges_read, "out.vtk", "all\n\nscattering:\n20 30\n70");
  std::system((GetPathToInstallDir() + "iges2vtk log.txt").c_str());  // NOLINT
  std::string log = GetFileContent("log.txt");
  ASSERT_THAT(log.find("log:\n"), Ne(std::string::npos));
  ASSERT_THAT(log.find(std::string("The splines at positions 0, 1 in file ") + iges_read
                           + " have been written to out.vtk"), Ne(std::string::npos));
  remove("log.txt");
  remove("out.vtk");
}

TEST_F(Iges2vtkExecutable, ConvertsFileCorrectly) {  // NOLINT
  CreateLogFile(iges_read, "out.vtk", "all\n\nscattering:\n20 30\n70");
  std::system((GetPathToInstallDir() + "iges2vtk log.txt").c_str());  // NOLINT
  std::string vtk_file = GetFileContent("out.vtk");
  ASSERT_THAT(vtk_file.find("# vtk DataFile Version 3.0\nSpline from Splinelib\nASCII\n"), Ne(std::string::npos));
  ASSERT_THAT(vtk_file.find("DATASET UNSTRUCTURED_GRID\nPOINTS 722 double\n"), Ne(std::string::npos));
  ASSERT_THAT(vtk_file.find("CELLS 670 3210\n"), Ne(std::string::npos));
  ASSERT_THAT(vtk_file.find("CELL_TYPES 670\n"), Ne(std::string::npos));
  remove("log.txt");
  remove("out.vtk");
}

class Iges2xmlExecutable : public Test {
 public:
  Iges2xmlExecutable() = default;

 protected:
  io::IGESReader iges_reader_;
  io::XMLReader xml_reader_;
};

TEST_F(Iges2xmlExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "iges2xml -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "iges2xml --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

TEST_F(Iges2xmlExecutable, WritesLog) {  // NOLINT
  CreateLogFile(iges_read_2, "out.xml", "1");
  std::system((GetPathToInstallDir() + "iges2xml log.txt").c_str());  // NOLINT
  std::string log = GetFileContent("log.txt");
  ASSERT_THAT(log.find("log:\n"), Ne(std::string::npos));
  ASSERT_THAT(log.find(std::string("The spline at position 1 in file ") + iges_read_2
                           + " has been written to out.xml"), Ne(std::string::npos));
  remove("log.txt");
  remove("out.xml");
}

TEST_F(Iges2xmlExecutable, ConvertsFileCorrectly) {  // NOLINT
  CreateLogFile(iges_read_2, "out.xml");
  std::vector<std::any> iges_splines = iges_reader_.ReadFile(iges_read_2);
  auto iges_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_splines[0]);
  auto iges_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(iges_splines[1]);
  std::system((GetPathToInstallDir() + "iges2xml log.txt").c_str());  // NOLINT
  std::vector<std::any> xml_splines = xml_reader_.ReadFile("out.xml");
  auto xml_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_splines[0]);
  auto xml_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(xml_splines[1]);
  ASSERT_THAT(xml_bspline_2d->AreEqual(*iges_bspline_2d.get(), 1e-6), true);
  ASSERT_THAT(xml_nurbs_1d->AreEqual(*iges_nurbs_1d.get(), 1e-6), true);
  remove("log.txt");
  remove("out.xml");
}

class Irit2igesExecutable : public Test {
 public:
  Irit2igesExecutable() = default;

 protected:
  io::IGESReader iges_reader_;
  io::IRITReader irit_reader_;
};

TEST_F(Irit2igesExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "irit2iges -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "irit2iges --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

TEST_F(Irit2igesExecutable, WritesLog) {  // NOLINT
  CreateLogFile(path_to_irit_file, "out.iges");
  std::system((GetPathToInstallDir() + "irit2iges log.txt").c_str());  // NOLINT
  std::string log = GetFileContent("log.txt");
  ASSERT_THAT(log.find("log:\n"), Ne(std::string::npos));
  ASSERT_THAT(log.find(std::string("The splines at positions 0, 1, 2, 3 in file ") + path_to_irit_file
                           + " have been written to out.iges"), Ne(std::string::npos));
  remove("log.txt");
  remove("out.iges");
}

TEST_F(Irit2igesExecutable, ConvertsFileCorrectly) {  // NOLINT
  CreateLogFile(path_to_irit_file, "out.iges", "1\n2\n3");
  std::vector<std::any> irit_splines = irit_reader_.ReadFile(path_to_irit_file);
  auto irit_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(irit_splines[1]);
  auto irit_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(irit_splines[2]);
  auto irit_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[3]);
  std::system((GetPathToInstallDir() + "irit2iges log.txt").c_str());  // NOLINT
  std::vector<std::any> iges_splines = iges_reader_.ReadFile("out.iges");
  auto iges_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(iges_splines[0]);
  auto iges_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_splines[1]);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[2]);
  ASSERT_THAT(irit_nurbs_1d->AreGeometricallyEqual(*iges_nurbs_1d.get()), true);
  ASSERT_THAT(irit_bspline_2d->AreEqual(*iges_bspline_2d.get()), true);
  ASSERT_THAT(irit_nurbs_2d->AreGeometricallyEqual(*iges_nurbs_2d.get()), true);
  remove("log.txt");
  remove("out.iges");
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

TEST_F(Irit2vtkExecutable, WritesLog) {  // NOLINT
  CreateLogFile(path_to_irit_file, "out.vtk", "0\n3\n5\n\nscattering:\n30\n10 20\n6 7 5");
  std::system((GetPathToInstallDir() + "irit2vtk log.txt").c_str());  // NOLINT
  std::string log = GetFileContent("log.txt");
  ASSERT_THAT(log.find("log:\n"), Ne(std::string::npos));
  ASSERT_THAT(log.find(std::string("The splines at positions 0, 3, 5 in file ") + path_to_irit_file
                           + " have been written to out.vtk"), Ne(std::string::npos));
  remove("log.txt");
  remove("out.vtk");
}

TEST_F(Irit2vtkExecutable, ConvertsFileCorrectly) {  // NOLINT
  CreateLogFile(path_to_irit_file, "out.vtk", "all\n\nscattering:\n30\n40\n10 5\n10 20\n5 8 5\n9 4 8");
  std::system((GetPathToInstallDir() + "irit2vtk log.txt").c_str());  // NOLINT
  std::string vtk_file = GetFileContent("out.vtk");
  ASSERT_THAT(vtk_file.find("# vtk DataFile Version 3.0\nSpline from Splinelib\nASCII\n"), Ne(std::string::npos));
  ASSERT_THAT(vtk_file.find("DATASET UNSTRUCTURED_GRID\nPOINTS 1143 double\n"), Ne(std::string::npos));
  ASSERT_THAT(vtk_file.find("CELLS 808 5852\n"), Ne(std::string::npos));
  ASSERT_THAT(vtk_file.find("CELL_TYPES 808\n"), Ne(std::string::npos));
  remove("log.txt");
  remove("out.vtk");
}

class Irit2xmlExecutable : public Test {
 public:
  Irit2xmlExecutable() = default;

 protected:
  io::IRITReader irit_reader_;
  io::XMLReader xml_reader_;
};

TEST_F(Irit2xmlExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "irit2xml -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "irit2xml --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

TEST_F(Irit2xmlExecutable, WritesLog) {  // NOLINT
  CreateLogFile(path_to_irit_file, "out.xml");
  std::system((GetPathToInstallDir() + "irit2xml log.txt").c_str());  // NOLINT
  std::string log = GetFileContent("log.txt");
  ASSERT_THAT(log.find("log:\n"), Ne(std::string::npos));
  ASSERT_THAT(log.find(std::string("The splines at positions 0, 1, 2, 3, 4, 5 in file ") + path_to_irit_file
                           + " have been written to out.xml"), Ne(std::string::npos));
  remove("log.txt");
  remove("out.xml");
}

TEST_F(Irit2xmlExecutable, ConvertsFileCorrectly) {  // NOLINT
  CreateLogFile(path_to_irit_file, "out.xml", "1\n2\n4");
  std::vector<std::any> irit_splines = irit_reader_.ReadFile(path_to_irit_file);
  auto irit_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(irit_splines[1]);
  auto irit_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(irit_splines[2]);
  auto irit_bspline_3d = std::any_cast<std::shared_ptr<spl::BSpline<3>>>(irit_splines[4]);
  std::system((GetPathToInstallDir() + "irit2xml log.txt").c_str());  // NOLINT
  std::vector<std::any> xml_splines = xml_reader_.ReadFile("out.xml");
  auto xml_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(xml_splines[0]);
  auto xml_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_splines[1]);
  auto xml_bspline_3d = std::any_cast<std::shared_ptr<spl::BSpline<3>>>(xml_splines[2]);
  ASSERT_THAT(irit_nurbs_1d->AreGeometricallyEqual(*xml_nurbs_1d.get()), true);
  ASSERT_THAT(irit_bspline_2d->AreEqual(*xml_bspline_2d.get()), true);
  ASSERT_THAT(irit_bspline_3d->AreEqual(*xml_bspline_3d.get()), true);
  remove("log.txt");
  remove("out.xml");
}

class Xml2igesExecutable : public Test {
 public:
  Xml2igesExecutable() = default;

 protected:
  io::IGESReader iges_reader_;
  io::XMLReader xml_reader_;
};

TEST_F(Xml2igesExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "xml2iges -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "xml2iges --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

TEST_F(Xml2igesExecutable, WritesLog) {  // NOLINT
  CreateLogFile(path_to_xml_file, "out.iges");
  std::system((GetPathToInstallDir() + "xml2iges log.txt").c_str());  // NOLINT
  std::string log = GetFileContent("log.txt");
  ASSERT_THAT(log.find("log:\n"), Ne(std::string::npos));
  ASSERT_THAT(log.find(std::string("The splines at positions 0, 1 in file ") + path_to_xml_file
                           + " have been written to out.iges"), Ne(std::string::npos));
  remove("log.txt");
  remove("out.iges");
}

TEST_F(Xml2igesExecutable, ConvertsFileCorrectly) {  // NOLINT
  CreateLogFile(path_to_xml_file, "out.iges");
  std::vector<std::any> xml_splines = xml_reader_.ReadFile(path_to_xml_file);
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines[0]);
  auto xml_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_splines[1]);
  std::system((GetPathToInstallDir() + "xml2iges log.txt").c_str());  // NOLINT
  std::vector<std::any> iges_splines = iges_reader_.ReadFile("out.iges");
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[0]);
  auto iges_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_splines[1]);
  ASSERT_THAT(xml_nurbs_2d->AreEqual(*iges_nurbs_2d.get()), true);
  ASSERT_THAT(xml_bspline_2d->AreEqual(*iges_bspline_2d.get()), true);
  remove("log.txt");
  remove("out.iges");
}

class Xml2iritExecutable : public Test {
 public:
  Xml2iritExecutable() = default;

 protected:
  io::IRITReader irit_reader_;
  io::XMLReader xml_reader_;
};

TEST_F(Xml2iritExecutable, PrintsHelp) {  // NOLINT
  std::string output = GetCommandOutput(GetPathToInstallDir() + "xml2irit -h");
  ASSERT_THAT(CompareToHelpOutput(output), true);
  output = GetCommandOutput(GetPathToInstallDir() + "xml2irit --help");
  ASSERT_THAT(CompareToHelpOutput(output), true);
}

TEST_F(Xml2iritExecutable, WritesLog) {  // NOLINT
  CreateLogFile(path_to_xml_file, "out.itd", "3");
  std::system((GetPathToInstallDir() + "xml2irit log.txt").c_str());  // NOLINT
  std::string log = GetFileContent("log.txt");
  ASSERT_THAT(log.find("log:\n"), Ne(std::string::npos));
  ASSERT_THAT(log.find(std::string("No splines in file ") + path_to_xml_file
                           + " have been written to out.itd"), Ne(std::string::npos));
  remove("log.txt");
  remove("out.itd");
}

TEST_F(Xml2iritExecutable, ConvertsFileCorrectly) {  // NOLINT
  CreateLogFile(path_to_xml_file, "out.itd");
  std::vector<std::any> xml_splines = xml_reader_.ReadFile(path_to_xml_file);
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines[0]);
  auto xml_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_splines[1]);
  std::system((GetPathToInstallDir() + "xml2irit log.txt").c_str());  // NOLINT
  std::vector<std::any> irit_splines = irit_reader_.ReadFile("out.itd");
  auto irit_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[0]);
  auto irit_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(irit_splines[1]);
  ASSERT_THAT(xml_nurbs_2d->AreEqual(*irit_nurbs_2d.get(), 1e-5), true);
  ASSERT_THAT(xml_bspline_2d->AreEqual(*irit_bspline_2d.get()), true);
  remove("log.txt");
  remove("out.itd");
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

TEST_F(Xml2vtkExecutable, WritesLog) {  // NOLINT
  CreateLogFile(path_to_xml_file, "out.vtk", "1\n3\n\nscattering:\n30 25\n3 4 5 6");
  std::system((GetPathToInstallDir() + "xml2vtk log.txt").c_str());  // NOLINT
  std::string log = GetFileContent("log.txt");
  ASSERT_THAT(log.find("log:\n"), Ne(std::string::npos));
  ASSERT_THAT(log.find(std::string("The spline at position 1 in file ") + path_to_xml_file
                           + " has been written to out.vtk"), Ne(std::string::npos));
  remove("log.txt");
  remove("out.vtk");
}

TEST_F(Xml2vtkExecutable, ConvertsFileCorrectly) {  // NOLINT
  CreateLogFile(path_to_xml_file, "out.vtk", "all\n\nscattering:\n10 5\n10 20");
  std::system((GetPathToInstallDir() + "xml2vtk log.txt").c_str());  // NOLINT
  std::string vtk_file = GetFileContent("out.vtk");
  ASSERT_THAT(vtk_file.find("# vtk DataFile Version 3.0\nSpline from Splinelib\nASCII\n"), Ne(std::string::npos));
  ASSERT_THAT(vtk_file.find("DATASET UNSTRUCTURED_GRID\nPOINTS 297 double\n"), Ne(std::string::npos));
  ASSERT_THAT(vtk_file.find("CELLS 250 1250\n"), Ne(std::string::npos));
  ASSERT_THAT(vtk_file.find("CELL_TYPES 250\n"), Ne(std::string::npos));
  remove("log.txt");
  remove("out.vtk");
}
