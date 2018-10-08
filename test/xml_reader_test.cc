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
#include <fstream>

#include "gmock/gmock.h"

#include "xml_reader.h"
#include "irit_reader.h"
#include "irit_writer.h"

using testing::Test;
using testing::DoubleEq;

class A2DSplineXMLReader : public Test {
 public:
  A2DSplineXMLReader() : xml_reader(std::make_unique<io::XMLReader<2>>()) {}

 protected:
  std::unique_ptr<io::XMLReader<2>> xml_reader;
};

TEST_F(A2DSplineXMLReader, ThrowsExceptionForNonExistingFile) {  // NOLINT
  ASSERT_THROW(xml_reader->ReadXMLFile("test.xml"), std::runtime_error);
}

TEST_F(A2DSplineXMLReader, GetsCorrectDegreeOfFirstSplineInFirstDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[0])->GetDegree(0).get(), 2);
}

TEST_F(A2DSplineXMLReader, GetsCorrectDegreeOfFirstSplineInSecondDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[0])->GetDegree(1).get(), 2);
}

TEST_F(A2DSplineXMLReader, GetsCorrectKnotOfFirstSplineInFirstDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[0])->GetKnotVector(0)->GetKnot(3).get(), DoubleEq(0.0625));
}

TEST_F(A2DSplineXMLReader, GetsCorrectKnotOfFirstSplineInSecondDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[0])->GetKnotVector(1)->GetKnot(3).get(), DoubleEq(0.125));
}

TEST_F(A2DSplineXMLReader, EvaluatesFirstSplineCorrectly) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[0])->Evaluate({ParamCoord(1), ParamCoord(1)}, {1})[0], DoubleEq(1));
}

TEST_F(A2DSplineXMLReader, GetsCorrectDegreeOfSecondSplineInFirstDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[1])->GetDegree(0).get(), 2);
}

TEST_F(A2DSplineXMLReader, GetsCorrectDegreeOfSecondSplineInSecondDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[1])->GetDegree(1).get(), 2);
}

TEST_F(A2DSplineXMLReader, GetsCorrectKnotOfSecondSplineInFirstDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[1])->GetKnotVector(0)->GetKnot(2).get(), DoubleEq(0.0));
}

TEST_F(A2DSplineXMLReader, GetsCorrectKnotOfSecondSplineInSecondDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[1])->GetKnotVector(1)->GetKnot(2).get(), DoubleEq(0.0));
}

TEST_F(A2DSplineXMLReader, EvaluatesSecondSplineCorrectly) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[1])->Evaluate({ParamCoord(0), ParamCoord(0)}, {0})[0], DoubleEq(-1));
}

TEST_F(A2DSplineXMLReader, ReturnsSameValuesBeforeAndAfterConvertingIRITToXMLFile) {  // NOLINT
  io::IRITWriter<2> irit_writer;
  irit_writer.ConvertXMLFileToIRITFile(path_to_xml_file, "converted_irit_file.itd");
  io::IRITReader<2> irit_reader;
  std::vector<std::any> spline_vector = irit_reader.ReadIRITFile("converted_irit_file.itd");
  ASSERT_THAT(spline_vector.size(), 2);

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      spline_vector[1])->Evaluate({ParamCoord(0), ParamCoord(0)}, {0})[0], DoubleEq(-1));

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      spline_vector[0])->Evaluate({ParamCoord(1), ParamCoord(1)}, {1})[0], DoubleEq(1));
  remove("converted_irit_file.itd");
}
