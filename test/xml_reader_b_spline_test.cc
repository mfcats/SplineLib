/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "xml_generator_b_spline.h"

#include <fstream>

#include "gmock/gmock.h"

using testing::Test;
using testing::DoubleEq;

class ABSplineXMLReader : public Test {
 public:
  ABSplineXMLReader() {
    xml_generator = std::make_unique<spl::XMLGenerator_B_Spline<2>>();
  }

 protected:
  const char *file = "spline_tank.xml";
  std::unique_ptr<spl::XMLGenerator_B_Spline<2>> xml_generator;
};

TEST_F(ABSplineXMLReader, OpensFile) {
  xml_generator->ReadXMLFile(file);
  ASSERT_THAT(xml_generator->ReadXMLFile(file)->GetDegree(0), 2);
  ASSERT_THAT(xml_generator->ReadXMLFile(file)->GetDegree(1), 2);
  ASSERT_THAT(xml_generator->ReadXMLFile(file)->GetKnotVector(0).GetKnot(3).get(), DoubleEq(0.0625));
  ASSERT_THAT(xml_generator->ReadXMLFile(file)->GetKnotVector(1).GetKnot(3).get(), DoubleEq(0.125));
  ASSERT_THAT(xml_generator->ReadXMLFile(file)->Evaluate({ParamCoord(1), ParamCoord(1)}, {1})[0], DoubleEq(1));
}
