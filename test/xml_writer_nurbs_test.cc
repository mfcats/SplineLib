/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "xml_writer_nurbs.h"

#include <fstream>

#include "gmock/gmock.h"

using testing::Test;

class ANURBSXMLWriter : public Test {
 public:
  ANURBSXMLWriter() {
    knot_vector = {std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5},
                                                                      ParamCoord{1}, ParamCoord{1}}))};
    degree_ = {Degree{1}};
    control_points_ = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0}))
    };
    std::vector<double> weights = {2.0, 1.75, 0.36};
    physical_space = std::make_shared<spl::WeightedPhysicalSpace<1>>(control_points_, weights, std::array<int, 1>({3}));
    parameter_space = std::make_shared<spl::ParameterSpace<1>>(knot_vector, degree_);
    std::vector<spl::NURBS<1>> splines;
    splines.emplace_back(parameter_space, physical_space);
    xml_writer = std::make_unique<io::XMLWriterNURBS<1>>(splines);
  }

 protected:
  std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector;
  std::unique_ptr<io::XMLWriterNURBS<1>> xml_writer;
  std::shared_ptr<spl::WeightedPhysicalSpace<1>> physical_space;
  std::shared_ptr<spl::ParameterSpace<1>> parameter_space;
  std::array<Degree, 1> degree_;
  std::vector<baf::ControlPoint> control_points_;
};

TEST_F(ANURBSXMLWriter, IsCreated) {  // NOLINT
  xml_writer->WriteXMLFile("nurbs.xml");
  std::ifstream newFile;
  newFile.open("nurbs.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("nurbs.xml");
}

TEST_F(ANURBSXMLWriter, CreatesCorrectXMLFile) {  // NOLINT
  xml_writer->WriteXMLFile("nurbs.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("nurbs.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("nurbs.xml");
}

TEST_F(ANURBSXMLWriter, CreatesSplineList) {  // NOLINT
  xml_writer->WriteXMLFile("nurbs.xml");
  pugi::xml_document doc;
  doc.load_file("nurbs.xml");
  ASSERT_STREQ(doc.first_child().name(), "SplineList");
  remove("nurbs.xml");
}

TEST_F(ANURBSXMLWriter, CreatesSpline) {  // NOLINT
  xml_writer->WriteXMLFile("nurbs.xml");
  pugi::xml_document doc;
  doc.load_file("nurbs.xml");
  ASSERT_STREQ(doc.first_child().first_child().name(), "SplineEntry");
  remove("nurbs.xml");
}

TEST_F(ANURBSXMLWriter, CreatesWeights) {  // NOLINT
  xml_writer->WriteXMLFile("nurbs.xml");
  pugi::xml_document doc;
  doc.load_file("nurbs.xml");
  ASSERT_STREQ(doc.first_child().first_child().child("wght").name(), "wght");
  remove("nurbs.xml");
}

class A2DNURBSXMLWriter : public Test {
 public:
  A2DNURBSXMLWriter() {
    knot_vector = {std::make_unique<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5},
                                                                      ParamCoord{1}, ParamCoord{1}})),
                   std::make_unique<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5},
                                                                      ParamCoord{1}, ParamCoord{1}}))};
    degree_ = {Degree{1}, Degree{2}};
    control_points_ = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0}))
    };
    std::vector<double> weights = {2.0, 1.75, 0.36, 1.0, 1.0, 0.05};
    physical_space =
        std::make_shared<spl::WeightedPhysicalSpace<2>>(control_points_, weights, std::array<int, 2>({3, 2}));
    parameter_space = std::make_shared<spl::ParameterSpace<2>>(knot_vector, degree_);
    std::vector<spl::NURBS<2>> splines;
    splines.emplace_back(parameter_space, physical_space);
    xml_writer = std::make_unique<io::XMLWriterNURBS<2>>(splines);
  }

 protected:
  std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector;
  std::unique_ptr<io::XMLWriterNURBS<2>> xml_writer;
  std::shared_ptr<spl::WeightedPhysicalSpace<2>> physical_space;
  std::shared_ptr<spl::ParameterSpace<2>> parameter_space;
  std::array<Degree, 2> degree_;
  std::vector<baf::ControlPoint> control_points_;
};

TEST_F(A2DNURBSXMLWriter, IsCreated) {  // NOLINT
  xml_writer->WriteXMLFile("2d_nurbs.xml");
  std::ifstream newFile;
  newFile.open("2d_nurbs.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("2d_nurbs.xml");
}

TEST_F(A2DNURBSXMLWriter, CreatesCorrectXMLFile) {  // NOLINT
  xml_writer->WriteXMLFile("2d_nurbs.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("2d_nurbs.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("2d_nurbs.xml");
}
