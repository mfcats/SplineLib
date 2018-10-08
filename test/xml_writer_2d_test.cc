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

class A2DBSplineForXML {  // NOLINT
 public:
  A2DBSplineForXML() {
    std::array<Degree, 2> degree = {Degree{1}, Degree{2}};
    std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector = {
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0, 2.0}))
    };
    b_spline_ = std::make_shared<spl::BSpline<2>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> b_spline_;
  virtual ~A2DBSplineForXML() = default;
};

class A2DNURBSForXML {  // NOLINT
 public:
  A2DNURBSForXML() {
    std::array<Degree, 2> degree = {Degree{1}, Degree{2}};
    std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector =
        {std::make_unique<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5},
                                                            ParamCoord{1}, ParamCoord{1}})),
         std::make_unique<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                            ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.6})),
        baf::ControlPoint(std::vector<double>({1.0, 0.8})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0}))
    };
    std::vector<double> weights = {2.0, 1.75, 0.36, 1.0, 1.0, 0.05, 1.0, 1.0, 1.0};

    nurbs_ = std::make_shared<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<2>> nurbs_;
  virtual ~A2DNURBSForXML() = default;
};

class A2DXMLWriter : public Test, public A2DBSplineForXML, public A2DNURBSForXML {
 public:
  A2DXMLWriter() : xml_writer_(std::make_unique<io::XMLWriter<2>>()) {
    std::any b_spline_any = std::make_any<std::shared_ptr<spl::BSpline<2>>>(b_spline_);
    std::any nurbs_any = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbs_);
    splines_ = {b_spline_any, nurbs_any};
  }

 protected:
  std::unique_ptr<io::XMLWriter<2>> xml_writer_;
  std::vector<std::any> splines_;
};

TEST_F(A2DXMLWriter, IsCreated) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "2d_splines.xml");
  std::ifstream newFile;
  newFile.open("2d_splines.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("2d_splines.xml");
}

TEST_F(A2DXMLWriter, CreatesCorrectXMLFile) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "2d_splines.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("2d_splines.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("2d_splines.xml");
}

TEST_F(A2DXMLWriter, CreatesSplineListWith2Entries) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "2d_splines.xml");
  pugi::xml_document doc;
  doc.load_file("2d_splines.xml");
  ASSERT_STREQ(doc.first_child().name(), "SplineList");
  ASSERT_STREQ(doc.first_child().attribute("NumberOfSplines").value(), "2");
  remove("2d_splines.xml");
}

TEST_F(A2DXMLWriter, Creates2SplineEntries) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "2d_splines.xml");
  pugi::xml_document doc;
  doc.load_file("2d_splines.xml");
  ASSERT_STREQ(doc.child("SplineList").first_child().name(), "SplineEntry");
  ASSERT_STREQ(doc.child("SplineList").first_child().next_sibling().name(), "SplineEntry");
  remove("2d_splines.xml");
}

TEST_F(A2DXMLWriter, CreatesNoWeightsForBSpline) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "2d_splines.xml");
  pugi::xml_document doc;
  doc.load_file("2d_splines.xml");
  ASSERT_STREQ(doc.child("SplineList").child("SplineEntry").child("wght").name(), "");
  ASSERT_STREQ(doc.child("SplineList").child("SplineEntry").next_sibling().child("wght").name(), "wght");
  remove("2d_splines.xml");
}

TEST_F(A2DXMLWriter, WritesSplineDimension2ForBothSplines) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "2d_splines.xml");
  pugi::xml_document doc;
  doc.load_file("2d_splines.xml");
  pugi::xml_node spline_node = doc.child("SplineList").child("SplineEntry");
  for (int i = 0; i < 2; i++, spline_node = spline_node.next_sibling()) {
    ASSERT_STREQ(spline_node.attribute("splDim").value(), "2");
  }
  remove("2d_splines.xml");
}

TEST_F(A2DXMLWriter, WritesCorrectSpaceDimensions) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "2d_splines.xml");
  pugi::xml_document doc;
  doc.load_file("2d_splines.xml");
  ASSERT_STREQ(doc.child("SplineList").child("SplineEntry").attribute("spaceDim").value(), "3");
  ASSERT_STREQ(doc.child("SplineList").child("SplineEntry").next_sibling().attribute("spaceDim").value(), "2");
  remove("2d_splines.xml");
}

TEST_F(A2DXMLWriter, ReturnsSameValuesBeforeAndAfterWritingAndReadingXMLFile) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "2d_splines.xml");
  std::unique_ptr<io::XMLReader<2>> xml_reader(std::make_unique<io::XMLReader<2>>());
  auto bspline_after = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_reader->ReadXMLFile("2d_splines.xml")[0]);
  auto nurbs_after = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_reader->ReadXMLFile("2d_splines.xml")[1]);

  ASSERT_THAT(b_spline_->Evaluate({ParamCoord(0.47681), ParamCoord(0.685409)}, {0})[0],
              DoubleEq(bspline_after->Evaluate({ParamCoord(0.47681), ParamCoord(0.685409)}, {0})[0]));
  ASSERT_THAT(b_spline_->Evaluate({ParamCoord(0.47681), ParamCoord(0.685409)}, {1})[0],
              DoubleEq(bspline_after->Evaluate({ParamCoord(0.47681), ParamCoord(0.685409)}, {1})[0]));

  ASSERT_THAT(nurbs_->Evaluate({ParamCoord(0.13697), ParamCoord(0.33246)}, {0})[0],
              DoubleEq(nurbs_after->Evaluate({ParamCoord(0.13697), ParamCoord(0.33246)}, {0})[0]));
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord(0.13697), ParamCoord(0.33246)}, {1})[0],
              DoubleEq(nurbs_after->Evaluate({ParamCoord(0.13697), ParamCoord(0.33246)}, {1})[0]));
  remove("2d_splines.xml");
}
