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

#include "xml_reader.h"

using testing::Test;
using testing::DoubleEq;

class ABSplineXMLWriter : public Test {
 public:
  ABSplineXMLWriter() {
    std::array<Degree, 1> degree = {Degree{1}};
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0}))
    };
    b_spline_ = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
    std::any b_spline_any = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_);
    std::vector<std::any> splines = {b_spline_any};
    xml_writer_ = std::make_unique<io::XMLWriter<1>>(splines);
  }

 protected:
  std::unique_ptr<io::XMLWriter<1>> xml_writer_;
  std::shared_ptr<spl::BSpline<1>> b_spline_;
};

TEST_F(ABSplineXMLWriter, IsCreated) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  std::ifstream newFile;
  newFile.open("bspline.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, CreatesCorrectXMLFile) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("bspline.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, CreatesSplineList) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.first_child().name(), "SplineList");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, CreatesSpline) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.first_child().first_child().name(), "SplineEntry");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, CreatesNoWeights) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.first_child().first_child().child("wght").name(), "");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, HasSplineDimension1) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.child("SplineList").child("SplineEntry").attribute("splDim").value(), "1");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, HasSpaceDimension2) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.child("SplineList").child("SplineEntry").attribute("spaceDim").value(), "2");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriter, ReturnsSameValuesBeforeAndAfterWritingAndReadingXMLFile) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  std::unique_ptr<io::XMLReader<1>> xml_reader(std::make_unique<io::XMLReader<1>>());
  auto bspline_after = std::any_cast<spl::BSpline<1>>(xml_reader->ReadXMLFile("bspline.xml")[0]);
  ASSERT_THAT(b_spline_->Evaluate({ParamCoord(0.75839)}, {0})[0],
              DoubleEq(bspline_after.Evaluate({ParamCoord(0.75839)}, {0})[0]));
}

class ABSplineXMLWriterWithSpaceDimension1 : public Test {
 public:
  ABSplineXMLWriterWithSpaceDimension1() {
    std::array<Degree, 1> degree = {Degree{1}};
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0})),
        baf::ControlPoint(std::vector<double>({0.5})),
        baf::ControlPoint(std::vector<double>({1.0}))
    };
    b_spline_ = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
    std::any b_spline_any = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_);
    std::vector<std::any> splines = {b_spline_any};
    xml_writer_ = std::make_unique<io::XMLWriter<1>>(splines);
  }

 protected:
  std::unique_ptr<io::XMLWriter<1>> xml_writer_;
  std::shared_ptr<spl::BSpline<1>> b_spline_;
};

TEST_F(ABSplineXMLWriterWithSpaceDimension1, IsCreated) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  std::ifstream newFile;
  newFile.open("bspline.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriterWithSpaceDimension1, CreatesCorrectXMLFile) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("bspline.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriterWithSpaceDimension1, HasSpaceDimension1) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.child("SplineList").child("SplineEntry").attribute("spaceDim").value(), "1");
  remove("bspline.xml");
}

class ABSplineXMLWriterWithSpaceDimension3 : public Test {
 public:
  ABSplineXMLWriterWithSpaceDimension3() {
    std::array<Degree, 1> degree = {Degree{1}};
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0}))
    };
    b_spline_ = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
    std::any b_spline_any = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_);
    std::vector<std::any> splines = {b_spline_any};
    xml_writer_ = std::make_unique<io::XMLWriter<1>>(splines);
  }

 protected:
  std::unique_ptr<io::XMLWriter<1>> xml_writer_;
  std::shared_ptr<spl::BSpline<1>> b_spline_;
};

TEST_F(ABSplineXMLWriterWithSpaceDimension3, IsCreated) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  std::ifstream newFile;
  newFile.open("bspline.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriterWithSpaceDimension3, CreatesCorrectXMLFile) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("bspline.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("bspline.xml");
}

TEST_F(ABSplineXMLWriterWithSpaceDimension3, HasSpaceDimension3) {  // NOLINT
  xml_writer_->WriteXMLFile("bspline.xml");
  pugi::xml_document doc;
  doc.load_file("bspline.xml");
  ASSERT_STREQ(doc.child("SplineList").child("SplineEntry").attribute("spaceDim").value(), "3");
  remove("bspline.xml");
}

class A2DBSplineXMLWriter : public Test {
 public:
  A2DBSplineXMLWriter() {
    std::array<Degree, 2> degree = {Degree{1}, Degree{2}};
    std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector = {
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0}))
    };
    b_spline_ = std::make_shared<spl::BSpline<2>>(knot_vector, degree, control_points);
    std::any b_spline_any = std::make_any<std::shared_ptr<spl::BSpline<2>>>(b_spline_);
    splines_ = {b_spline_any};
    xml_writer_ = std::make_unique<io::XMLWriter<2>>(splines_);
  }

 protected:
  std::unique_ptr<io::XMLWriter<2>> xml_writer_;
  std::shared_ptr<spl::BSpline<2>> b_spline_;
  std::vector<std::any> splines_;
};

TEST_F(A2DBSplineXMLWriter, IsCreated) {  // NOLINT
  xml_writer_->WriteXMLFile("2d_bspline.xml");
  std::ifstream newFile;
  newFile.open("2d_bspline.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("2d_bspline.xml");
}

TEST_F(A2DBSplineXMLWriter, CreatesCorrectXMLFile) {  // NOLINT
  xml_writer_->WriteXMLFile("2d_bspline.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("2d_bspline.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("2d_bspline.xml");
}

TEST_F(A2DBSplineXMLWriter, HasSplineDimension2) {  // NOLINT
  xml_writer_->WriteXMLFile("2d_bspline.xml");
  pugi::xml_document doc;
  doc.load_file("2d_bspline.xml");
  ASSERT_STREQ(doc.child("SplineList").child("SplineEntry").attribute("splDim").value(), "2");
  remove("2d_bspline.xml");
}

TEST_F(A2DBSplineXMLWriter, CanWriteTwoSplines) {  // NOLINT
  splines_.emplace_back(b_spline_);
  xml_writer_ = std::make_unique<io::XMLWriter<2>>(splines_);
  xml_writer_->WriteXMLFile("2d_bspline.xml");
  pugi::xml_document doc;
  doc.load_file("2d_bspline.xml");
  ASSERT_STREQ(doc.child("SplineList").attribute("NumberOfSplines").value(), "2");
  remove("2d_bspline.xml");
}
