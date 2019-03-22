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

#include "gmock/gmock.h"

#include "io_converter.h"

using testing::Test;
using testing::DoubleNear;
using testing::Ne;

class AnIOConverter : public Test {
 public:
  AnIOConverter() : io_converter_(std::make_unique<io::IOConverter>()) {}

 protected:
  std::unique_ptr<io::IOConverter> io_converter_;
  io::IGESReader iges_reader_;
  io::IRITReader irit_reader_;
  io::XMLReader xml_reader_;
};

TEST_F(AnIOConverter, ReturnsSameValueBeforeAndAfterConvertingSplinesFromIGESFileToIRITFile) {  // NOLINT
  std::vector<std::any> iges_splines = iges_reader_.ReadFile(iges_read);
  ASSERT_THAT(iges_splines.size(), 2);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[0]);
  auto iges_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines[1]);
  io_converter_->ConvertFile(iges_read, "converted_irit_file.itd");
  std::vector<std::any> irit_splines = irit_reader_.ReadFile("converted_irit_file.itd");
  ASSERT_THAT(irit_splines.size(), iges_splines.size());
  auto irit_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[0]);
  auto irit_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(irit_splines[1]);
  ASSERT_THAT(irit_nurbs_2d->AreEqual(*iges_nurbs_2d.get(), 1e-6), true);
  ASSERT_THAT(irit_bspline_1d->AreEqual(*iges_bspline_1d.get(), 1e-6), true);
  remove("converted_irit_file.itd");
}

TEST_F(AnIOConverter, ReturnsSameValueBeforeAndAfterConvertingSplinesFromIGESFileToIRITFile2) {  // NOLINT
  std::vector<std::any> iges_splines = iges_reader_.ReadFile(iges_read_2);
  ASSERT_THAT(iges_splines.size(), 2);
  auto iges_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_splines[0]);
  auto iges_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(iges_splines[1]);
  io_converter_->ConvertFile(iges_read_2, "converted_irit_file.itd");
  std::vector<std::any> irit_splines = irit_reader_.ReadFile("converted_irit_file.itd");
  ASSERT_THAT(irit_splines.size(), iges_splines.size());
  auto irit_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(irit_splines[0]);
  auto irit_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(irit_splines[1]);
  ASSERT_THAT(irit_bspline_2d->AreEqual(*iges_bspline_2d.get()), true);
  ASSERT_THAT(irit_nurbs_1d->AreEqual(*iges_nurbs_1d.get(), 1e-6), true);
  remove("converted_irit_file.itd");
}

TEST_F(AnIOConverter, ReturnsSameValueBeforeAndAfterConvertingSplinesFromIGESFileToXMLFile) {  // NOLINT
  std::vector<std::any> iges_splines = iges_reader_.ReadFile(iges_read);
  ASSERT_THAT(iges_splines.size(), 2);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[0]);
  auto iges_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines[1]);
  io_converter_->ConvertFile(iges_read, "converted_xml_file.xml");
  std::vector<std::any> xml_splines = xml_reader_.ReadFile("converted_xml_file.xml");
  ASSERT_THAT(xml_splines.size(), iges_splines.size());
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines[0]);
  auto xml_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(xml_splines[1]);
  ASSERT_THAT(xml_nurbs_2d->AreEqual(*iges_nurbs_2d.get(), 1e-6), true);
  ASSERT_THAT(xml_bspline_1d->AreEqual(*iges_bspline_1d.get(), 1e-6), true);
  remove("converted_xml_file.xml");
}

TEST_F(AnIOConverter, ConvertsSplinesFromIGESFileToVTKFile) {  // NOLINT
  std::vector<std::any> iges_splines = iges_reader_.ReadFile(iges_read);
  ASSERT_THAT(iges_splines.size(), 2);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[0]);
  auto iges_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines[1]);
  io_converter_->ConvertFile(iges_read, "converted_vtk_file.vtk", {{20, 30}, {70}});
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

TEST_F(AnIOConverter, ReturnsSameValueBeforeAndAfterConvertingSplinesFromXMLFileToIGESFile) {  // NOLINT
  std::vector<std::any> xml_splines = xml_reader_.ReadFile(path_to_xml_file);
  ASSERT_THAT(xml_splines.size(), 4);
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines[0]);
  auto xml_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_splines[1]);
  io_converter_->ConvertFile(path_to_xml_file, "converted_iges_file.iges");
  std::vector<std::any> iges_splines = iges_reader_.ReadFile("converted_iges_file.iges");
  ASSERT_THAT(iges_splines.size(), xml_splines.size() - 2);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[0]);
  auto iges_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_splines[1]);
  ASSERT_THAT(xml_nurbs_2d->AreEqual(*iges_nurbs_2d.get()), true);
  ASSERT_THAT(xml_bspline_2d->AreEqual(*iges_bspline_2d.get()), true);
  remove("converted_iges_file.iges");
}

TEST_F(AnIOConverter, ReturnsSameValueBeforeAndAfterConvertingSplinesFromXMLFileToIRITFile) {  // NOLINT
  std::vector<std::any> xml_splines = xml_reader_.ReadFile(path_to_xml_file);
  ASSERT_THAT(xml_splines.size(), 4);
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines[0]);
  auto xml_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_splines[1]);
  io_converter_->ConvertFile(path_to_xml_file, "converted_irit_file.itd");
  std::vector<std::any> irit_splines = irit_reader_.ReadFile("converted_irit_file.itd");
  ASSERT_THAT(irit_splines.size(), xml_splines.size() - 2);
  auto irit_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[0]);
  auto irit_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(irit_splines[1]);
  ASSERT_THAT(xml_nurbs_2d->AreEqual(*irit_nurbs_2d.get(), 1e-5), true);
  ASSERT_THAT(xml_bspline_2d->AreEqual(*irit_bspline_2d.get()), true);
  remove("converted_irit_file.itd");
}

TEST_F(AnIOConverter, ReturnsSameValueBeforeAndAfterConvertingSplinesFromIRITFileToIGESFile) {  // NOLINT
  std::vector<std::any> irit_splines = irit_reader_.ReadFile(path_to_irit_file);
  ASSERT_THAT(irit_splines.size(), 6);
  auto irit_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(irit_splines[0]);
  auto irit_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(irit_splines[1]);
  auto irit_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(irit_splines[2]);
  auto irit_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[3]);
  io_converter_->ConvertFile(path_to_irit_file, "converted_iges_file.iges");
  std::vector<std::any> iges_splines = iges_reader_.ReadFile("converted_iges_file.iges");
  ASSERT_THAT(iges_splines.size(), irit_splines.size() - 2);
  auto iges_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines[0]);
  auto iges_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(iges_splines[1]);
  auto iges_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_splines[2]);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[3]);
  ASSERT_THAT(irit_bspline_1d->AreEqual(*iges_bspline_1d.get()), true);
  ASSERT_THAT(irit_nurbs_1d->AreGeometricallyEqual(*iges_nurbs_1d.get()), true);
  ASSERT_THAT(irit_bspline_2d->AreEqual(*iges_bspline_2d.get()), true);
  ASSERT_THAT(irit_nurbs_2d->AreGeometricallyEqual(*iges_nurbs_2d.get()), true);
  remove("converted_iges_file.iges");
}

TEST_F(AnIOConverter, ReturnsSameValueBeforeAndAfterConvertingSplinesFromIRITFileToXMLFile) {  // NOLINT
  std::vector<std::any> irit_splines = irit_reader_.ReadFile(path_to_irit_file);
  ASSERT_THAT(irit_splines.size(), 6);
  auto irit_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(irit_splines[0]);
  auto irit_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(irit_splines[1]);
  auto irit_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(irit_splines[2]);
  auto irit_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[3]);
  auto irit_bspline_3d = std::any_cast<std::shared_ptr<spl::BSpline<3>>>(irit_splines[4]);
  auto irit_nurbs_3d = std::any_cast<std::shared_ptr<spl::NURBS<3>>>(irit_splines[5]);
  io_converter_->ConvertFile(path_to_irit_file, "converted_xml_file.xml");
  std::vector<std::any> xml_splines = xml_reader_.ReadFile("converted_xml_file.xml");
  ASSERT_THAT(xml_splines.size(), irit_splines.size());
  auto xml_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(xml_splines[0]);
  auto xml_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(xml_splines[1]);
  auto xml_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_splines[2]);
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines[3]);
  auto xml_bspline_3d = std::any_cast<std::shared_ptr<spl::BSpline<3>>>(xml_splines[4]);
  auto xml_nurbs_3d = std::any_cast<std::shared_ptr<spl::NURBS<3>>>(xml_splines[5]);
  ASSERT_THAT(xml_bspline_1d->AreEqual(*irit_bspline_1d.get(), 1e-6), true);
  ASSERT_THAT(xml_nurbs_1d->AreEqual(*irit_nurbs_1d.get()), true);
  ASSERT_THAT(xml_bspline_2d->AreEqual(*irit_bspline_2d.get()), true);
  ASSERT_THAT(xml_nurbs_2d->AreEqual(*irit_nurbs_2d.get()), true);
  ASSERT_THAT(xml_bspline_3d->AreEqual(*irit_bspline_3d.get()), true);
  ASSERT_THAT(xml_nurbs_3d->AreEqual(*irit_nurbs_3d.get()), true);
  remove("converted_xml_file.xml");
}

TEST_F(AnIOConverter, ThrowsForWrongTypeOfInputFile) {  // NOLINT
  ASSERT_THROW(io_converter_->ConvertFile("file.txt", "file.iges"), std::runtime_error);
}

TEST_F(AnIOConverter, ThrowsForWrongTypeOfOutputFile) {  // NOLINT
  ASSERT_THROW(io_converter_->ConvertFile(path_to_xml_file, "file.txt"), std::runtime_error);
}
