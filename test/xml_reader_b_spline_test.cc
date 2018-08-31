/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "xml_reader.h"

#include <config.h>
#include <fstream>

#include "gmock/gmock.h"

using testing::Test;
using testing::DoubleEq;

class ABSplineXMLReader : public Test {
 public:
  ABSplineXMLReader() : xml_reader(std::make_unique<io::XMLReader<2>>()) {}

 protected:
  std::unique_ptr<io::XMLReader<2>> xml_reader;
};

TEST_F(ABSplineXMLReader, OpensFile) {  // NOLINT
  xml_reader->ReadXMLFile(path_to_xml_file);
  ASSERT_THAT(xml_reader->ReadXMLFile(path_to_xml_file)->GetDegree(0), 2);
  ASSERT_THAT(xml_reader->ReadXMLFile(path_to_xml_file)->GetDegree(1), 2);
  ASSERT_THAT(xml_reader->ReadXMLFile(path_to_xml_file)->GetKnotVector(0).GetKnot(3).get(), DoubleEq(0.0625));
  ASSERT_THAT(xml_reader->ReadXMLFile(path_to_xml_file)->GetKnotVector(1).GetKnot(3).get(), DoubleEq(0.125));
  ASSERT_THAT(xml_reader->ReadXMLFile(path_to_xml_file)->Evaluate({ParamCoord(1), ParamCoord(1)}, {1})[0], DoubleEq(1));
}
