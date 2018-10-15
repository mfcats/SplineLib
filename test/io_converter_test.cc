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
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[0]);
  auto iges_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines[1]);
  io_converter_->ConvertFile(iges_read, "converted_irit_file.itd");
  std::vector<std::any> irit_splines = irit_reader_.ReadFile("converted_irit_file.itd");
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[0]);
  auto xml_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(irit_splines[1]);
  ASSERT_THAT(xml_nurbs_2d->Evaluate({ParamCoord(0.34867)}, {0})[0],
              DoubleNear(iges_nurbs_2d->Evaluate({ParamCoord(0.34867)}, {0})[0], 0.00001));
  ASSERT_THAT(xml_nurbs_2d->Evaluate({ParamCoord(0.34867)}, {1})[0],
              DoubleNear(iges_nurbs_2d->Evaluate({ParamCoord(0.34867)}, {1})[0], 0.00001));

  ASSERT_THAT(xml_bspline_1d->Evaluate({ParamCoord(0.34867)}, {0})[0],
              DoubleNear(iges_bspline_1d->Evaluate({ParamCoord(0.34867)}, {0})[0], 0.00001));
  remove("converted_irit_file.itd");
}

TEST_F(AnIOConverter, ReturnsSameValueBeforeAndAfterConvertingSplinesFromIGESFileToXMLFile) {  // NOLINT
  std::vector<std::any> iges_splines = iges_reader_.ReadFile(iges_read);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[0]);
  auto iges_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines[1]);
  io_converter_->ConvertFile(iges_read, "converted_xml_file.xml");
  std::vector<std::any> xml_splines = xml_reader_.ReadFile("converted_xml_file.xml");
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines[0]);
  auto xml_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(xml_splines[1]);
  ASSERT_THAT(xml_nurbs_2d->Evaluate({ParamCoord(0.76584)}, {0})[0],
              DoubleNear(iges_nurbs_2d->Evaluate({ParamCoord(0.76584)}, {0})[0], 0.00001));
  ASSERT_THAT(xml_nurbs_2d->Evaluate({ParamCoord(0.76584)}, {1})[0],
              DoubleNear(iges_nurbs_2d->Evaluate({ParamCoord(0.76584)}, {1})[0], 0.00001));

  ASSERT_THAT(xml_bspline_1d->Evaluate({ParamCoord(0.76584)}, {0})[0],
              DoubleNear(iges_bspline_1d->Evaluate({ParamCoord(0.76584)}, {0})[0], 0.00001));
  remove("converted_xml_file.xml");
}

TEST_F(AnIOConverter, ReturnsSameValueBeforeAndAfterConvertingSplinesFromXMLFileToIGESFile) {  // NOLINT
  std::vector<std::any> xml_splines = xml_reader_.ReadFile(path_to_xml_file);
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines[0]);
  auto xml_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_splines[1]);
  io_converter_->ConvertFile(path_to_xml_file, "converted_iges_file.iges");
  std::vector<std::any> iges_splines = iges_reader_.ReadFile("converted_iges_file.iges");
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[0]);
  auto iges_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_splines[1]);
  ASSERT_THAT(xml_nurbs_2d->Evaluate({ParamCoord(0.00124)}, {0})[0],
              DoubleNear(iges_nurbs_2d->Evaluate({ParamCoord(0.00124)}, {0})[0], 0.00001));
  ASSERT_THAT(xml_nurbs_2d->Evaluate({ParamCoord(0.00124)}, {1})[0],
              DoubleNear(iges_nurbs_2d->Evaluate({ParamCoord(0.00124)}, {1})[0], 0.00001));

  ASSERT_THAT(xml_bspline_2d->Evaluate({ParamCoord(0.00124)}, {0})[0],
              DoubleNear(iges_bspline_2d->Evaluate({ParamCoord(0.00124)}, {0})[0], 0.00001));
  ASSERT_THAT(xml_bspline_2d->Evaluate({ParamCoord(0.00124)}, {1})[0],
              DoubleNear(iges_bspline_2d->Evaluate({ParamCoord(0.00124)}, {1})[0], 0.00001));
  remove("converted_iges_file.iges");
}

TEST_F(AnIOConverter, ReturnsSameValueBeforeAndAfterConvertingSplinesFromXMLFileToIRITFile) {  // NOLINT
  std::vector<std::any> xml_splines = xml_reader_.ReadFile(path_to_xml_file);
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines[0]);
  auto xml_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_splines[1]);
  io_converter_->ConvertFile(path_to_xml_file, "converted_irit_file.itd");
  std::vector<std::any> irit_splines = irit_reader_.ReadFile("converted_irit_file.itd");
  ASSERT_THAT(irit_splines.size(), 2);
  auto irit_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[0]);
  auto irit_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(irit_splines[1]);
  ASSERT_THAT(xml_nurbs_2d->Evaluate({ParamCoord(0.99979)}, {0})[0],
              DoubleNear(irit_nurbs_2d->Evaluate({ParamCoord(0.99979)}, {0})[0], 0.00001));
  ASSERT_THAT(xml_nurbs_2d->Evaluate({ParamCoord(0.99979)}, {1})[0],
              DoubleNear(irit_nurbs_2d->Evaluate({ParamCoord(0.99979)}, {1})[0], 0.00001));

  ASSERT_THAT(xml_bspline_2d->Evaluate({ParamCoord(0.99979)}, {0})[0],
              DoubleNear(irit_bspline_2d->Evaluate({ParamCoord(0.99979)}, {0})[0], 0.00001));
  ASSERT_THAT(xml_bspline_2d->Evaluate({ParamCoord(0.99979)}, {1})[0],
              DoubleNear(irit_bspline_2d->Evaluate({ParamCoord(0.99979)}, {1})[0], 0.00001));
  remove("converted_irit_file.itd");
}

TEST_F(AnIOConverter, ThrowsWhenConvertingSplinesFromIRITFileToIGESFile) {  // NOLINT
  std::vector<std::any> irit_splines = irit_reader_.ReadFile(path_to_iris_file);
  auto irit_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(irit_splines[0]);
  auto irit_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(irit_splines[1]);
  auto irit_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(irit_splines[2]);
  auto irit_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[3]);
  io_converter_->ConvertFile(path_to_iris_file, "converted_iges_file.iges");
  std::vector<std::any> iges_splines = iges_reader_.ReadFile("converted_iges_file.iges");
  auto iges_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines[0]);
  auto iges_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(iges_splines[1]);
  auto iges_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_splines[2]);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines[3]);
  ASSERT_THAT(iges_splines.size(), 4);
  ASSERT_THAT(irit_bspline_1d->Evaluate({ParamCoord(0.00124)}, {0})[0],
              DoubleNear(iges_bspline_1d->Evaluate({ParamCoord(0.00124)}, {0})[0], 0.00001));
  ASSERT_THAT(irit_bspline_1d->Evaluate({ParamCoord(0.00124)}, {1})[0],
              DoubleNear(iges_bspline_1d->Evaluate({ParamCoord(0.00124)}, {1})[0], 0.00001));
  ASSERT_THAT(irit_nurbs_2d->Evaluate({ParamCoord(0.00124)}, {0})[0],
              DoubleNear(iges_nurbs_2d->Evaluate({ParamCoord(0.00124)}, {0})[0], 0.00001));
  ASSERT_THAT(irit_nurbs_2d->Evaluate({ParamCoord(0.00124)}, {1})[0],
              DoubleNear(iges_nurbs_2d->Evaluate({ParamCoord(0.00124)}, {1})[0], 0.00001));
  remove("converted_iges_file.iges");
}

TEST_F(AnIOConverter, ReturnsSameValueBeforeAndAfterConvertingSplinesFromIRITFileToXMLFile) {  // NOLINT
  std::vector<std::any> irit_splines = irit_reader_.ReadFile(path_to_iris_file);
  auto irit_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(irit_splines[0]);
  auto irit_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_splines[3]);
  auto irit_nurbs_3d = std::any_cast<std::shared_ptr<spl::NURBS<3>>>(irit_splines[5]);
  io_converter_->ConvertFile(path_to_iris_file, "converted_xml_file.xml");
  std::vector<std::any> xml_splines = xml_reader_.ReadFile("converted_xml_file.xml");
  auto xml_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(xml_splines[0]);
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines[3]);
  auto xml_nurbs_3d = std::any_cast<std::shared_ptr<spl::NURBS<3>>>(xml_splines[5]);
  ASSERT_THAT(xml_bspline_1d->Evaluate({ParamCoord(0.76584)}, {0})[0],
              DoubleNear(irit_bspline_1d->Evaluate({ParamCoord(0.76584)}, {0})[0], 0.00001));

  ASSERT_THAT(xml_nurbs_2d->Evaluate({ParamCoord(0.76584)}, {0})[0],
              DoubleNear(irit_nurbs_2d->Evaluate({ParamCoord(0.76584)}, {0})[0], 0.00001));
  ASSERT_THAT(xml_nurbs_2d->Evaluate({ParamCoord(0.76584)}, {1})[0],
              DoubleNear(irit_nurbs_2d->Evaluate({ParamCoord(0.76584)}, {1})[0], 0.00001));

  ASSERT_THAT(xml_nurbs_3d->Evaluate({ParamCoord(0.76584)}, {0})[0],
              DoubleNear(irit_nurbs_3d->Evaluate({ParamCoord(0.76584)}, {0})[0], 0.00001));
  ASSERT_THAT(xml_nurbs_3d->Evaluate({ParamCoord(0.76584)}, {1})[0],
              DoubleNear(irit_nurbs_3d->Evaluate({ParamCoord(0.76584)}, {1})[0], 0.00001));
  ASSERT_THAT(xml_nurbs_3d->Evaluate({ParamCoord(0.76584)}, {2})[0],
              DoubleNear(irit_nurbs_3d->Evaluate({ParamCoord(0.76584)}, {2})[0], 0.00001));
  remove("converted_xml_file.xml");
}
