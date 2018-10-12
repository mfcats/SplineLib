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
  AnIOConverter() : io_converter(std::make_unique<io::IOConverter>()) {}

 protected:
  std::unique_ptr<io::IOConverter> io_converter;
};

TEST_F(AnIOConverter, ReturnsSameValueBeforeAndAfterConvertingSplinesFromIGESFileToXMLFile) {  // NOLINT
  io::IGESReader iges_reader;
  std::vector<std::any> iges_splines_1 = iges_reader.ReadFile(iges_read);
  std::vector<std::any> iges_splines_2 = iges_reader.ReadFile(iges_read_2);
  auto iges_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_splines_1[0]);
  auto iges_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_splines_1[1]);
  auto iges_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_splines_2[0]);
  auto iges_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(iges_splines_2[1]);
  io_converter->ConvertIGESFileToXMLFile(iges_read, "converted_xml_file_1.xml");
  io_converter->ConvertIGESFileToXMLFile(iges_read_2, "converted_xml_file_2.xml");
  io::XMLReader xml_reader;
  std::vector<std::any> xml_splines_1 = xml_reader.ReadFile("converted_xml_file_1.xml");
  std::vector<std::any> xml_splines_2 = xml_reader.ReadFile("converted_xml_file_2.xml");
  auto xml_nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_splines_1[0]);
  auto xml_bspline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(xml_splines_1[1]);
  auto xml_bspline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_splines_2[0]);
  auto xml_nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(xml_splines_2[1]);
  ASSERT_THAT(xml_nurbs_2d->Evaluate({ParamCoord(0.76584)}, {1})[0],
              DoubleNear(iges_nurbs_2d->Evaluate({ParamCoord(0.76584)}, {1})[0], 0.00001));
  ASSERT_THAT(xml_nurbs_2d->Evaluate({ParamCoord(0.76584)}, {1})[0],
              DoubleNear(iges_nurbs_2d->Evaluate({ParamCoord(0.76584)}, {1})[0], 0.00001));

  ASSERT_THAT(xml_bspline_1d->Evaluate({ParamCoord(0.76584)}, {0})[0],
              DoubleNear(iges_bspline_1d->Evaluate({ParamCoord(0.76584)}, {0})[0], 0.00001));

  ASSERT_THAT(xml_bspline_2d->Evaluate({ParamCoord(0.76584)}, {1})[0],
              DoubleNear(iges_bspline_2d->Evaluate({ParamCoord(0.76584)}, {1})[0], 0.00001));
  ASSERT_THAT(xml_bspline_2d->Evaluate({ParamCoord(0.76584)}, {1})[0],
              DoubleNear(iges_bspline_2d->Evaluate({ParamCoord(0.76584)}, {1})[0], 0.00001));

  ASSERT_THAT(xml_nurbs_1d->Evaluate({ParamCoord(0.76584)}, {0})[0],
              DoubleNear(iges_nurbs_1d->Evaluate({ParamCoord(0.76584)}, {0})[0], 0.00001));
  remove("converted_xml_file_1.xml");
  remove("converted_xml_file_2.xml");
}
