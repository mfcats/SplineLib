/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "xml_writer.h"

#include <fstream>

#include "gmock/gmock.h"

using testing::Test;

class ANURBSXMLWriter : public Test {
 public:
  ANURBSXMLWriter() {
    std::array<Degree, 1> degree = {Degree{1}};
    std::array<std::shared_ptr<baf::KnotVector>, 1>
        knot_vector = {std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5},
                                                                          ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0}))
    };
    std::vector<double> weights = {2.0, 1.75, 0.36};
    nurbs_ = std::make_shared<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
    std::any b_spline_any = std::make_any<std::shared_ptr<spl::NURBS<1>>>(nurbs_);
    std::vector<std::any> splines = {b_spline_any};
    xml_writer_ = std::make_unique<io::XMLWriter<1>>(splines);
  }

 protected:
  std::unique_ptr<io::XMLWriter<1>> xml_writer_;
  std::shared_ptr<spl::NURBS<1>> nurbs_;
};

TEST_F(ANURBSXMLWriter, IsCreated) {  // NOLINT
  xml_writer_->WriteXMLFile("nurbs.xml");
  std::ifstream newFile;
  newFile.open("nurbs.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("nurbs.xml");
}

TEST_F(ANURBSXMLWriter, CreatesCorrectXMLFile) {  // NOLINT
  xml_writer_->WriteXMLFile("nurbs.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("nurbs.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("nurbs.xml");
}

TEST_F(ANURBSXMLWriter, CreatesSplineList) {  // NOLINT
  xml_writer_->WriteXMLFile("nurbs.xml");
  pugi::xml_document doc;
  doc.load_file("nurbs.xml");
  ASSERT_STREQ(doc.first_child().name(), "SplineList");
  remove("nurbs.xml");
}

TEST_F(ANURBSXMLWriter, CreatesSpline) {  // NOLINT
  xml_writer_->WriteXMLFile("nurbs.xml");
  pugi::xml_document doc;
  doc.load_file("nurbs.xml");
  ASSERT_STREQ(doc.first_child().first_child().name(), "SplineEntry");
  remove("nurbs.xml");
}

TEST_F(ANURBSXMLWriter, CreatesWeights) {  // NOLINT
  xml_writer_->WriteXMLFile("nurbs.xml");
  pugi::xml_document doc;
  doc.load_file("nurbs.xml");
  ASSERT_STREQ(doc.first_child().first_child().child("wght").name(), "wght");
  remove("nurbs.xml");
}

class A2DNURBSXMLWriter : public Test {
 public:
  A2DNURBSXMLWriter() {
    std::array<Degree, 2> degree = {Degree{1}, Degree{2}};
    std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector =
        {std::make_unique<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5},
                                                            ParamCoord{1}, ParamCoord{1}})),
         std::make_unique<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5},
                                                            ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0}))
    };
    std::vector<double> weights = {2.0, 1.75, 0.36, 1.0, 1.0, 0.05};
    nurbs_ = std::make_shared<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
    std::any nurbs_any = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbs_);
    std::vector<std::any> splines = {nurbs_any};
    xml_writer_ = std::make_unique<io::XMLWriter<2>>(splines);
  }

 protected:
  std::unique_ptr<io::XMLWriter<2>> xml_writer_;
  std::shared_ptr<spl::NURBS<2>> nurbs_;
};

TEST_F(A2DNURBSXMLWriter, IsCreated) {  // NOLINT
  xml_writer_->WriteXMLFile("2d_nurbs.xml");
  std::ifstream newFile;
  newFile.open("2d_nurbs.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("2d_nurbs.xml");
}

TEST_F(A2DNURBSXMLWriter, CreatesCorrectXMLFile) {  // NOLINT
  xml_writer_->WriteXMLFile("2d_nurbs.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("2d_nurbs.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("2d_nurbs.xml");
}
