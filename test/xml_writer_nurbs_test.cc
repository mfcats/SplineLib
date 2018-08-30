/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <fstream>

#include "gmock/gmock.h"

#include "pugixml.hpp"

#include "xml_generator_nurbs.h"

using testing::Test;

class ANURBSXMLWriter : public Test {
 public:
  ANURBSXMLWriter() {
    std::array<baf::KnotVector, 1> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}})};
    degree_ = {1};
    control_points_ = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0}))
    };
    std::vector<double> weights = {2.0, 1.75, 0.36};
    physical_space = spl::WeightedPhysicalSpace<1>(control_points_, weights, {3});
    parameter_space = spl::ParameterSpace<1>(knot_vector, degree_);
    xml_generator = std::make_unique<spl::XMLGenerator_NURBS<1>>(physical_space, parameter_space);
  }

 protected:
  std::unique_ptr<spl::XMLGenerator_NURBS<1>> xml_generator;
  spl::WeightedPhysicalSpace<1> physical_space;
  spl::ParameterSpace<1> parameter_space;
  std::array<int, 1> degree_;
  std::vector<baf::ControlPoint> control_points_;
};

TEST_F(ANURBSXMLWriter, IsCreated) {
  xml_generator->WriteXMLFile("nurbs.xml");
  std::ifstream newFile;
  newFile.open("nurbs.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
}

TEST_F(ANURBSXMLWriter, CreatesCorrectXMLFile) {
  xml_generator->WriteXMLFile("nurbs.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("nurbs.xml");
  ASSERT_STREQ(result.description(), "No error");
}

TEST_F(ANURBSXMLWriter, CreatesSplineList) {
  xml_generator->WriteXMLFile("nurbs.xml");
  pugi::xml_document doc;
  doc.load_file("nurbs.xml");
  ASSERT_STREQ(doc.first_child().name(), "SplineList");
}

TEST_F(ANURBSXMLWriter, CreatesSpline) {
  xml_generator->WriteXMLFile("nurbs.xml");
  pugi::xml_document doc;
  doc.load_file("nurbs.xml");
  ASSERT_STREQ(doc.first_child().first_child().name(), "SplineEntry");
}
