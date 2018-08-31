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
using testing::Eq;

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

TEST_F(ABSplineXMLWriter, IsCreated) {  // NOLINT
  xml_generator->WriteXMLFile("bspline.xml");
  std::ifstream newFile;
  newFile.open("bspline.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, CreatesCorrectXMLFile) {  // NOLINT
  xml_generator->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("bspline.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, CreatesSplineList) {  // NOLINT
  xml_generator->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.first_child().name(), "SplineList");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, CreatesSpline) {  // NOLINT
  xml_generator->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.first_child().first_child().name(), "SplineEntry");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, CreatesNoWeights) {  // NOLINT
  xml_generator->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.first_child().first_child().child("wght").name(), "");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, HasSplineDimension1) {  // NOLINT
  xml_generator->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.child("SplineList").child("SplineEntry").attribute("splDim").value(), "1");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, HasSpaceDimension2) {  // NOLINT
  xml_generator->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.child("SplineList").child("SplineEntry").attribute("spaceDim").value(), "2");
  remove("bspline.xml");
}

class ABSplineXMLWriterWithSpaceDimension3 : public Test {
 public:
  ABSplineXMLWriterWithSpaceDimension3() {
    std::array<baf::KnotVector, 1> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}})};
    degree_ = {1};
    control_points_ = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0}))
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

TEST_F(ABSplineXMLWriterWithSpaceDimension3, IsCreated) {  // NOLINT
  xml_generator->WriteXMLFile("bspline.xml");
  std::ifstream newFile;
  newFile.open("bspline.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriterWithSpaceDimension3, CreatesCorrectXMLFile) {  // NOLINT
  xml_generator->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("bspline.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriterWithSpaceDimension3, HasSpaceDimension3) {  // NOLINT
  xml_generator->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.child("SplineList").child("SplineEntry").attribute("spaceDim").value(), "3");
  remove("bspline.xml");
}

class A2DBSplineXMLWriter : public Test {
 public:
  A2DBSplineXMLWriter() {
    std::array<baf::KnotVector, 2> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}}),
         baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}})};
    degree_ = {1, 2};
    control_points_ = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0}))
    };
    physical_space = spl::PhysicalSpace<2>(control_points_, {3, 2});
    parameter_space = spl::ParameterSpace<2>(knot_vector, degree_);
    xml_generator = std::make_unique<spl::XMLGenerator_B_Spline<2>>(physical_space, parameter_space);
  }

 protected:
  std::unique_ptr<spl::XMLGenerator_B_Spline<2>> xml_generator;
  spl::PhysicalSpace<2> physical_space;
  spl::ParameterSpace<2> parameter_space;
  std::array<int, 2> degree_;
  std::vector<baf::ControlPoint> control_points_;
};

TEST_F(A2DBSplineXMLWriter, IsCreated) {  // NOLINT
  xml_generator->WriteXMLFile("2d_bspline.xml");
  std::ifstream newFile;
  newFile.open("2d_bspline.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("2d_bspline.xml");
}

TEST_F(A2DBSplineXMLWriter, CreatesCorrectXMLFile) {  // NOLINT
  xml_generator->WriteXMLFile("2d_bspline.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("2d_bspline.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("2d_bspline.xml");
}

TEST_F(A2DBSplineXMLWriter, HasSplineDimension2) {  // NOLINT
  xml_generator->WriteXMLFile("2d_bspline.xml");
  pugi::xml_document doc;
  doc.load_file("2d_bspline.xml");
  ASSERT_STREQ(doc.child("SplineList").child("SplineEntry").attribute("splDim").value(), "2");
  remove("2d_bspline.xml");
}