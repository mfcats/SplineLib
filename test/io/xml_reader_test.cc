/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <config_xml.h>
#include <fstream>

#include "gmock/gmock.h"

#include "xml_reader.h"
#include "irit_reader.h"
#include "irit_writer.h"

using testing::Test;
using testing::DoubleEq;

class AnXMLReader : public Test {
 public:
  AnXMLReader() : xml_reader(std::make_unique<io::XMLReader>()) {}

 protected:
  std::unique_ptr<io::XMLReader> xml_reader;
};

TEST_F(AnXMLReader, ThrowsExceptionForNonExistingFile) {  // NOLINT
  ASSERT_THROW(xml_reader->ReadFile("test.xml"), std::runtime_error);
}

TEST_F(AnXMLReader, FindsFourSplines) {  // NOLINT
  ASSERT_THAT(xml_reader->ReadFile(path_to_xml_file).size(), 4);
}

TEST_F(AnXMLReader, GetsCorrectDegrees) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[0])->GetDegree(0).get(), 2);
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[0])->GetDegree(1).get(), 2);

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[1])->GetDegree(0).get(), 2);
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[1])->GetDegree(1).get(), 2);

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[2])->GetDegree(0).get(), 2);
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[2])->GetDegree(1).get(), 1);
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[2])->GetDegree(2).get(), 1);
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[2])->GetDegree(3).get(), 1);

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[3])->GetDegree(0).get(), 1);
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[3])->GetDegree(1).get(), 1);
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[3])->GetDegree(2).get(), 1);
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[3])->GetDegree(3).get(), 1);
}

TEST_F(AnXMLReader, GetsCorrectKnots) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[0])->GetKnotVector(0)->GetKnot(3).get(), DoubleEq(0.0625));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[0])->GetKnotVector(1)->GetKnot(3).get(), DoubleEq(0.125));

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[1])->GetKnotVector(0)->GetKnot(2).get(), DoubleEq(0.0));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[1])->GetKnotVector(1)->GetKnot(2).get(), DoubleEq(0.0));

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[2])->GetKnotVector(0)->GetKnot(2).get(), DoubleEq(0.0));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[2])->GetKnotVector(1)->GetKnot(1).get(), DoubleEq(0.0));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[2])->GetKnotVector(2)->GetKnot(2).get(), DoubleEq(1.0));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[2])->GetKnotVector(3)->GetKnot(3).get(), DoubleEq(1.0));

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[3])->GetKnotVector(0)->GetKnot(2).get(), DoubleEq(1.0));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[3])->GetKnotVector(1)->GetKnot(1).get(), DoubleEq(0.0));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[3])->GetKnotVector(2)->GetKnot(2).get(), DoubleEq(1.0));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<4>>>(
      xml_reader->ReadFile(path_to_xml_file)[3])->GetKnotVector(3)->GetKnot(3).get(), DoubleEq(1.0));
}

TEST_F(AnXMLReader, EvaluatesSplinesCorrectly) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[0])->Evaluate({ParamCoord(1), ParamCoord(1)}, {0})[0], DoubleEq(0));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[0])->Evaluate({ParamCoord(1), ParamCoord(1)}, {1})[0], DoubleEq(1));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[0])->Evaluate({ParamCoord(1), ParamCoord(1)}, {2})[0], DoubleEq(4));

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[1])->Evaluate({ParamCoord(0), ParamCoord(0)}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[1])->Evaluate({ParamCoord(0), ParamCoord(0)}, {1})[0], DoubleEq(-2));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      xml_reader->ReadFile(path_to_xml_file)[1])->Evaluate({ParamCoord(0), ParamCoord(0)}, {2})[0], DoubleEq(3));

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<4>>>(xml_reader->ReadFile(path_to_xml_file)[2])->
      Evaluate({ParamCoord(0), ParamCoord(0), ParamCoord(0), ParamCoord(0)}, {0})[0], DoubleEq(-1));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<4>>>(xml_reader->ReadFile(path_to_xml_file)[2])->
      Evaluate({ParamCoord(0), ParamCoord(0), ParamCoord(0), ParamCoord(0)}, {1})[0], DoubleEq(-1));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<4>>>(xml_reader->ReadFile(path_to_xml_file)[2])->
      Evaluate({ParamCoord(0), ParamCoord(0), ParamCoord(0), ParamCoord(0)}, {2})[0], DoubleEq(0));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<4>>>(xml_reader->ReadFile(path_to_xml_file)[2])->
      Evaluate({ParamCoord(0), ParamCoord(0), ParamCoord(0), ParamCoord(0)}, {3})[0], DoubleEq(1));

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<4>>>(xml_reader->ReadFile(path_to_xml_file)[3])->
      Evaluate({ParamCoord(0), ParamCoord(0), ParamCoord(0), ParamCoord(0)}, {0})[0], DoubleEq(1));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<4>>>(xml_reader->ReadFile(path_to_xml_file)[3])->
      Evaluate({ParamCoord(0), ParamCoord(0), ParamCoord(0), ParamCoord(0)}, {1})[0], DoubleEq(0.5));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<4>>>(xml_reader->ReadFile(path_to_xml_file)[3])->
      Evaluate({ParamCoord(0), ParamCoord(0), ParamCoord(0), ParamCoord(0)}, {2})[0], DoubleEq(0.2));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<4>>>(xml_reader->ReadFile(path_to_xml_file)[3])->
      Evaluate({ParamCoord(0), ParamCoord(0), ParamCoord(0), ParamCoord(0)}, {3})[0], DoubleEq(0.8));
}
