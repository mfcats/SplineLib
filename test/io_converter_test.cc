/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <config.h>

#include "gmock/gmock.h"

#include "io_converter.h"
#include "xml_reader.h"

using testing::Test;
using testing::DoubleNear;

class A1DIOConverter : public Test {
 public:
  A1DIOConverter() : io_converter(std::make_unique<io::IOConverter<1>>()) {}

 protected:
  std::unique_ptr<io::IOConverter<1>> io_converter;
};

TEST_F(A1DIOConverter, ReturnsSameValueBeforeAndAfterConverting1DBSplineFromIGESFileToXMLFile) {  // NOLINT
  io::IGESReader iges_reader;
  std::vector<std::any> iges_splines = iges_reader.ReadIGESFile(iges_read);
  auto iges_spline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines[1]);

  io_converter->ConvertIGESFileToXMLFile(iges_read, "converted_xml_file_1d.xml");
  io::XMLReader<1> xml_reader_1d;
  std::vector<std::any> xml_splines = xml_reader_1d.ReadXMLFile("converted_xml_file_1d.xml");
  ASSERT_THAT(xml_splines.size(), 1);
  auto xml_spline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(xml_splines[0]);
  ASSERT_THAT(xml_spline_1d->Evaluate({ParamCoord(0.76584)}, {0})[0],
              DoubleNear(iges_spline_1d->Evaluate({ParamCoord(0.76584)}, {0})[0], 0.00001));

  remove("converted_xml_file_1d.xml");
}

TEST_F(A1DIOConverter, ReturnsSameValueBeforeAndAfterConverting1DNURBSFromIGESFileToXMLFile) {  // NOLINT
  io::IGESReader iges_reader;
  std::vector<std::any> iges_splines = iges_reader.ReadIGESFile(iges_read_2);
  auto iges_spline_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(iges_splines[1]);

  io_converter->ConvertIGESFileToXMLFile(iges_read_2, "converted_xml_file_1d.xml");
  io::XMLReader<1> xml_reader_1d;
  std::vector<std::any> xml_splines = xml_reader_1d.ReadXMLFile("converted_xml_file_1d.xml");
  ASSERT_THAT(xml_splines.size(), 1);
  auto xml_spline_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(xml_splines[0]);
  ASSERT_THAT(xml_spline_1d->Evaluate({ParamCoord(0.76584)}, {0})[0],
              DoubleNear(iges_spline_1d->Evaluate({ParamCoord(0.76584)}, {0})[0], 0.00001));

  remove("converted_xml_file_1d.xml");
}

class A2DIOConverter : public Test {
 public:
  A2DIOConverter() : io_converter(std::make_unique<io::IOConverter<2>>()) {}

 protected:
  std::unique_ptr<io::IOConverter<2>> io_converter;
};

TEST_F(A2DIOConverter, ReturnsSameValueBeforeAndAfterConverting2DBSplineFromIGESFileToXMLFile) {  // NOLINT
  io::IGESReader iges_reader;
  std::vector<std::any> iges_splines = iges_reader.ReadIGESFile(iges_read_2);
  auto iges_spline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_splines[0]);

  io_converter->ConvertIGESFileToXMLFile(iges_read_2, "converted_xml_file_2d.xml");
  io::XMLReader<2> xml_reader_2d;
  std::vector<std::any> xml_splines = xml_reader_2d.ReadXMLFile("converted_xml_file_2d.xml");
  ASSERT_THAT(xml_splines.size(), 1);
  auto xml_spline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_splines[0]);
  ASSERT_THAT(xml_spline_2d->Evaluate({ParamCoord(0.76584)}, {1})[0],
              DoubleNear(iges_spline_2d->Evaluate({ParamCoord(0.76584)}, {1})[0], 0.00001));
  ASSERT_THAT(xml_spline_2d->Evaluate({ParamCoord(0.76584)}, {1})[0],
              DoubleNear(iges_spline_2d->Evaluate({ParamCoord(0.76584)}, {1})[0], 0.00001));

  remove("converted_xml_file_2d.xml");
}

TEST_F(A2DIOConverter, ReturnsSameValueBeforeAndAfterConverting2DNURBSFromIGESFileToXMLFile) {  // NOLINT
  io::IGESReader iges_reader;
  std::vector<std::any> iges_splines = iges_reader.ReadIGESFile(iges_read);
  auto iges_spline_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[0]);

  io_converter->ConvertIGESFileToXMLFile(iges_read, "converted_xml_file_2d.xml");
  io::XMLReader<2> xml_reader_2d;
  std::vector<std::any> xml_splines = xml_reader_2d.ReadXMLFile("converted_xml_file_2d.xml");
  ASSERT_THAT(xml_splines.size(), 1);
  auto xml_spline_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines[0]);
  ASSERT_THAT(xml_spline_2d->Evaluate({ParamCoord(0.76584)}, {1})[0],
              DoubleNear(iges_spline_2d->Evaluate({ParamCoord(0.76584)}, {1})[0], 0.00001));
  ASSERT_THAT(xml_spline_2d->Evaluate({ParamCoord(0.76584)}, {1})[0],
              DoubleNear(iges_spline_2d->Evaluate({ParamCoord(0.76584)}, {1})[0], 0.00001));

  remove("converted_xml_file_2d.xml");
}
