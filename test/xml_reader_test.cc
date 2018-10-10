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

class A2DXMLReader : public Test {
 public:
  A2DXMLReader() : xml_reader(std::make_unique<io::XMLReader>()) {}

 protected:
  std::unique_ptr<io::XMLReader> xml_reader;
};

TEST_F(A2DXMLReader, ThrowsExceptionForNonExistingFile) {  // NOLINT
  ASSERT_THROW(xml_reader->ReadXMLFile("test.xml"), std::runtime_error);
}

TEST_F(A2DXMLReader, FindsTwoSplines) {  // NOLINT
  ASSERT_THAT(xml_reader->ReadXMLFile(path_to_xml_file).size(), 2);
}

TEST_F(A2DXMLReader, GetsCorrectDegreeOfFirstSplineInFirstDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[0])->GetDegree(0).get(), 2);
}

TEST_F(A2DXMLReader, GetsCorrectDegreeOfFirstSplineInSecondDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[0])->GetDegree(1).get(), 2);
}

TEST_F(A2DXMLReader, GetsCorrectKnotOfFirstSplineInFirstDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[0])->GetKnotVector(0)->GetKnot(3).get(), DoubleEq(0.0625));
}

TEST_F(A2DXMLReader, GetsCorrectKnotOfFirstSplineInSecondDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[0])->GetKnotVector(1)->GetKnot(3).get(), DoubleEq(0.125));
}

TEST_F(A2DXMLReader, EvaluatesFirstSplineCorrectly) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[0])->Evaluate({ParamCoord(1), ParamCoord(1)}, {1})[0], DoubleEq(1));
}

TEST_F(A2DXMLReader, GetsCorrectDegreeOfSecondSplineInFirstDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[1])->GetDegree(0).get(), 2);
}

TEST_F(A2DXMLReader, GetsCorrectDegreeOfSecondSplineInSecondDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[1])->GetDegree(1).get(), 2);
}

TEST_F(A2DXMLReader, GetsCorrectKnotOfSecondSplineInFirstDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[1])->GetKnotVector(0)->GetKnot(2).get(), DoubleEq(0.0));
}

TEST_F(A2DXMLReader, GetsCorrectKnotOfSecondSplineInSecondDirection) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[1])->GetKnotVector(1)->GetKnot(2).get(), DoubleEq(0.0));
}

TEST_F(A2DXMLReader, EvaluatesSecondSplineCorrectly) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadXMLFile(path_to_xml_file)[1])->Evaluate({ParamCoord(0), ParamCoord(0)}, {0})[0], DoubleEq(-1));
}
