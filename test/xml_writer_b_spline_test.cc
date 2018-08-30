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
#include <stdio.h>

#include "gmock/gmock.h"

#include "pugixml.hpp"

#include "xml_generator_b_spline.h"

using testing::Test;

class ABSplineXMLWriter : public Test {
 public:
  ABSplineXMLWriter() {
    std::array<baf::KnotVector, 1> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}})};
    degree_ = {1};
    control_points_ = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0}))
    };
    physical_space = spl::PhysicalSpace<1>(control_points_, {3});
    parameter_space = spl::ParameterSpace<1>(knot_vector, degree_);
    xml_generator = std::make_unique<spl::XMLGenerator_B_Spline<1>>(physical_space, parameter_space);
  }

 protected:
  std::unique_ptr<spl::XMLGenerator_B_Spline<1>> xml_generator;
  spl::PhysicalSpace<1> physical_space;
  spl::ParameterSpace<1> parameter_space;
  std::array<int, 1> degree_;
  std::vector<baf::ControlPoint> control_points_;
};

TEST_F(ABSplineXMLWriter, IsCreated) {
  xml_generator->WriteXMLFile("bspline.xml");
  std::ifstream newFile;
  newFile.open("bspline.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, CreatesCorrectXMLFile) {
  xml_generator->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("bspline.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, CreatesSplineList) {
  xml_generator->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.first_child().name(), "SplineList");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, CreatesSpline) {
  xml_generator->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.first_child().first_child().name(), "SplineEntry");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, CreatesNoWeights) {
  xml_generator->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.first_child().first_child().child("wght").name(), "");
  remove("bspline.xml");
}
