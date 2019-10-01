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

#include <fstream>

#include "gmock/gmock.h"

#include "io_converter.h"

using testing::Test;
using testing::Ne;

using namespace splinelib::src;

class AnIGEStoIRITConverter : public Test {
 public:
  AnIGEStoIRITConverter() : io_converter_(std::make_unique<io::IOConverter>(iges_read_2, "converted_irit_file.itd")) {}

 protected:
  std::unique_ptr<io::IOConverter> io_converter_;
  io::IGESReader iges_reader_;
  io::IRITReader irit_reader_;
};

TEST_F(AnIGEStoIRITConverter, CorrectlyConvertsSplines) {  // NOLINT
  std::vector<std::any> iges_splines = iges_reader_.ReadFile(iges_read_2);
  ASSERT_THAT(iges_splines.size(), 2);
  auto iges_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_splines[0]);
  auto iges_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(iges_splines[1]);
  ASSERT_THAT(io_converter_->ConvertFile().size(), 2);
  std::vector<std::any> irit_splines = irit_reader_.ReadFile("converted_irit_file.itd");
  ASSERT_THAT(irit_splines.size(), iges_splines.size());
  auto irit_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(irit_splines[0]);
  auto irit_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(irit_splines[1]);
  ASSERT_THAT(irit_bspline_2d->AreEqual(*iges_bspline_2d.get()), true);
  ASSERT_THAT(irit_nurbs_1d->AreEqual(*iges_nurbs_1d.get(), 1e-6), true);
  remove("converted_irit_file.itd");
}

class AnIGEStoVTKConverter : public Test {
 public:
  AnIGEStoVTKConverter() : io_converter_(std::make_unique<io::IOConverter>(iges_read, "converted_vtk_file.vtk")) {}

 protected:
  std::unique_ptr<io::IOConverter> io_converter_;
  io::IGESReader iges_reader_;
};

TEST_F(AnIGEStoVTKConverter, CorrectlyConvertsSplines) {  // NOLINT
  std::vector<std::any> iges_splines = iges_reader_.ReadFile(iges_read);
  ASSERT_THAT(iges_splines.size(), 2);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[0]);
  auto iges_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines[1]);
  ASSERT_THAT(io_converter_->ConvertFile({}, {{20, 30}, {70}}).size(), 2);
  std::ifstream newFile;
  newFile.open("converted_vtk_file.vtk");
  ASSERT_THAT(newFile.good(), true);
  std::string line, file;
  while (getline(newFile, line)) {
    file += line + "\n";
  }
  ASSERT_THAT(file.find("# vtk DataFile Version 3.0\nSpline from Splinelib\nASCII\n"), Ne(std::string::npos));
  ASSERT_THAT(file.find("DATASET UNSTRUCTURED_GRID\nPOINTS 722 double\n"), Ne(std::string::npos));
  ASSERT_THAT(file.find("CELLS 670 3210\n"), Ne(std::string::npos));
  ASSERT_THAT(file.find("CELL_TYPES 670\n"), Ne(std::string::npos));
  remove("converted_vtk_file.vtk");
}

class AnIGEStoXMLConverter : public Test {
 public:
  AnIGEStoXMLConverter() : io_converter_(std::make_unique<io::IOConverter>(iges_read, "converted_xml_file.xml")) {}

 protected:
  std::unique_ptr<io::IOConverter> io_converter_;
  io::IGESReader iges_reader_;
  io::XMLReader xml_reader_;
};

TEST_F(AnIGEStoXMLConverter, CorrectlyConvertsSplines) {  // NOLINT
  std::vector<std::any> iges_splines = iges_reader_.ReadFile(iges_read);
  ASSERT_THAT(iges_splines.size(), 2);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[0]);
  auto iges_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines[1]);
  ASSERT_THAT(io_converter_->ConvertFile().size(), 2);
  std::vector<std::any> xml_splines = xml_reader_.ReadFile("converted_xml_file.xml");
  ASSERT_THAT(xml_splines.size(), iges_splines.size());
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines[0]);
  auto xml_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(xml_splines[1]);
  ASSERT_THAT(xml_nurbs_2d->AreEqual(*iges_nurbs_2d.get(), 1e-6), true);
  ASSERT_THAT(xml_bspline_1d->AreEqual(*iges_bspline_1d.get(), 1e-6), true);
  remove("converted_xml_file.xml");
}

class AnIRITtoIGESConverter : public Test {
 public:
  AnIRITtoIGESConverter() : io_converter_(std::make_unique<io::IOConverter>(path_to_irit_file,
                                                                            "converted_iges_file.iges")) {}

 protected:
  std::unique_ptr<io::IOConverter> io_converter_;
  io::IGESReader iges_reader_;
  io::IRITReader irit_reader_;
};

TEST_F(AnIRITtoIGESConverter, CorrectlyConvertsSplines) {  // NOLINT
  std::vector<std::any> irit_splines = irit_reader_.ReadFile(path_to_irit_file);
  ASSERT_THAT(irit_splines.size(), 6);
  auto irit_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(irit_splines[0]);
  auto irit_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[3]);
  testing::internal::CaptureStderr();
  ASSERT_THAT(io_converter_->ConvertFile({0, 3}).size(), 2);
  std::vector<std::any> iges_splines = iges_reader_.ReadFile("converted_iges_file.iges");
  ASSERT_THAT(iges_splines.size(), 2);
  auto iges_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines[0]);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[1]);
  ASSERT_THAT(irit_bspline_1d->AreEqual(*iges_bspline_1d.get()), true);
  ASSERT_THAT(irit_nurbs_2d->AreGeometricallyEqual(*iges_nurbs_2d.get()), true);
  remove("converted_iges_file.iges");
}

class AnIRITtoVTKConverter : public Test {
 public:
  AnIRITtoVTKConverter() : io_converter_(std::make_unique<io::IOConverter>(path_to_irit_file,
                                                                           "converted_vtk_file.vtk")) {}

 protected:
  std::unique_ptr<io::IOConverter> io_converter_;
};

TEST_F(AnIRITtoVTKConverter, CorrectlyConvertsSplines) {  // NOLINT
  ASSERT_THAT(io_converter_->ConvertFile({5}, {{12, 13, 9}}).size(), 1);
  std::ifstream newFile;
  newFile.open("converted_vtk_file.vtk");
  ASSERT_THAT(newFile.good(), true);
  std::string line, file;
  while (getline(newFile, line)) {
    file += line + "\n";
  }
  ASSERT_THAT(file.find("# vtk DataFile Version 3.0\nSpline from Splinelib\nASCII\n"), Ne(std::string::npos));
  ASSERT_THAT(file.find("DATASET UNSTRUCTURED_GRID\nPOINTS 1820 double\n"), Ne(std::string::npos));
  ASSERT_THAT(file.find("CELLS 1404 12636\n"), Ne(std::string::npos));
  ASSERT_THAT(file.find("CELL_TYPES 1404\n"), Ne(std::string::npos));
  remove("converted_vtk_file.vtk");
}

class AnIRITtoXMLConverter : public Test {
 public:
  AnIRITtoXMLConverter() : io_converter_(std::make_unique<io::IOConverter>(path_to_irit_file,
                                                                           "converted_xml_file.xml")) {}

 protected:
  std::unique_ptr<io::IOConverter> io_converter_;
  io::IRITReader irit_reader_;
  io::XMLReader xml_reader_;
};

TEST_F(AnIRITtoXMLConverter, CorrectlyConvertsSplines) {  // NOLINT
  std::vector<std::any> irit_splines = irit_reader_.ReadFile(path_to_irit_file);
  ASSERT_THAT(irit_splines.size(), 6);
  auto irit_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(irit_splines[1]);
  auto irit_bspline_3d = std::any_cast<std::shared_ptr<spl::BSpline<3>>>(irit_splines[4]);
  ASSERT_THAT(io_converter_->ConvertFile({1, 4}).size(), 2);
  std::vector<std::any> xml_splines = xml_reader_.ReadFile("converted_xml_file.xml");
  ASSERT_THAT(xml_splines.size(), 2);
  auto xml_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(xml_splines[0]);
  auto xml_bspline_3d = std::any_cast<std::shared_ptr<spl::BSpline<3>>>(xml_splines[1]);
  ASSERT_THAT(xml_nurbs_1d->AreEqual(*irit_nurbs_1d.get()), true);
  ASSERT_THAT(xml_bspline_3d->AreEqual(*irit_bspline_3d.get()), true);
  remove("converted_xml_file.xml");
}

class AnXMLtoIGESConverter : public Test {
 public:
  AnXMLtoIGESConverter() : io_converter_(std::make_unique<io::IOConverter>(path_to_xml_file,
                                                                           "converted_iges_file.iges")) {}

 protected:
  std::unique_ptr<io::IOConverter> io_converter_;
  io::IGESReader iges_reader_;
  io::XMLReader xml_reader_;
};

TEST_F(AnXMLtoIGESConverter, CorrectlyConvertsSplines) {  // NOLINT
  std::vector<std::any> xml_splines = xml_reader_.ReadFile(path_to_xml_file);
  ASSERT_THAT(xml_splines.size(), 4);
  auto xml_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_splines[1]);
  ASSERT_THAT(io_converter_->ConvertFile({1, 2}).size(), 1);
  std::vector<std::any> iges_splines = iges_reader_.ReadFile("converted_iges_file.iges");
  ASSERT_THAT(iges_splines.size(), 1);
  auto iges_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_splines[0]);
  ASSERT_THAT(xml_bspline_2d->AreEqual(*iges_bspline_2d.get()), true);
  remove("converted_iges_file.iges");
}

class AnXMLtoIRITConverter : public Test {
 public:
  AnXMLtoIRITConverter() : io_converter_(std::make_unique<io::IOConverter>(path_to_xml_file,
                                                                           "converted_irit_file.itd")) {}

 protected:
  std::unique_ptr<io::IOConverter> io_converter_;
  io::IRITReader irit_reader_;
  io::XMLReader xml_reader_;
};

TEST_F(AnXMLtoIRITConverter, CorrectlyConvertsSplines) {  // NOLINT
  std::vector<std::any> xml_splines = xml_reader_.ReadFile(path_to_xml_file);
  ASSERT_THAT(xml_splines.size(), 4);
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines[0]);
  ASSERT_THAT(io_converter_->ConvertFile({0, 2, 3}).size(), 1);
  std::vector<std::any> irit_splines = irit_reader_.ReadFile("converted_irit_file.itd");
  ASSERT_THAT(irit_splines.size(), 1);
  auto irit_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[0]);
  ASSERT_THAT(xml_nurbs_2d->AreEqual(*irit_nurbs_2d.get(), 1e-5), true);
  remove("converted_irit_file.itd");
}

class AnIOConverterWithWrongInputFileFormat : public Test {
 public:
  AnIOConverterWithWrongInputFileFormat() : io_converter_(std::make_unique<io::IOConverter>("file.vtk", "file.iges")) {}

 protected:
  std::unique_ptr<io::IOConverter> io_converter_;
};

TEST_F(AnIOConverterWithWrongInputFileFormat, ThrowsWhenTryingToConvert) {  // NOLINT
  ASSERT_THROW(io_converter_->ConvertFile(), std::runtime_error);
}

class AnIOConverterWithWrongOutputFileFormat : public Test {
 public:
  AnIOConverterWithWrongOutputFileFormat() : io_converter_(std::make_unique<io::IOConverter>(iges_read, "file.txt")) {}

 protected:
  std::unique_ptr<io::IOConverter> io_converter_;
};

TEST_F(AnIOConverterWithWrongOutputFileFormat, ThrowsWhenTryingToConvert) {  // NOLINT
  ASSERT_THROW(io_converter_->ConvertFile(), std::runtime_error);
}
