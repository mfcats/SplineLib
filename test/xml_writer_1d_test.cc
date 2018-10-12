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
#include "config.in.h"

using testing::Test;
using testing::DoubleEq;

class A1DBSplineForXMLWithSpaceDim1 {  // NOLINT
 public:
  A1DBSplineForXMLWithSpaceDim1() {
    std::array<Degree, 1> degree = {Degree{1}};
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0})),
        baf::ControlPoint(std::vector<double>({0.5})),
        baf::ControlPoint(std::vector<double>({1.0}))
    };
    b_spline_1_ = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> b_spline_1_;
  virtual ~A1DBSplineForXMLWithSpaceDim1() = default;
};

class A1DBSplineForXMLWithSpaceDim2 {  // NOLINT
 public:
  A1DBSplineForXMLWithSpaceDim2() {
    std::array<Degree, 1> degree = {Degree{1}};
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0}))
    };
    b_spline_2_ = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> b_spline_2_;
  virtual ~A1DBSplineForXMLWithSpaceDim2() = default;
};

class A1DBSplineForXMLWithSpaceDim3 {  // NOLINT
 public:
  A1DBSplineForXMLWithSpaceDim3() {
    std::array<Degree, 1> degree = {Degree{1}};
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0}))
    };
    b_spline_3_ = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> b_spline_3_;
  virtual ~A1DBSplineForXMLWithSpaceDim3() = default;
};

class A1DNURBSForXML {  // NOLINT
 public:
  A1DNURBSForXML() {
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
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_;
  virtual ~A1DNURBSForXML() = default;
};

class A1DXMLWriter : public Test, public A1DBSplineForXMLWithSpaceDim1, public A1DBSplineForXMLWithSpaceDim2,
                     public A1DBSplineForXMLWithSpaceDim3, public A1DNURBSForXML {
 public:
  A1DXMLWriter() : xml_writer_(std::make_unique<io::XMLWriter>()) {
    std::any b_spline_1_any = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_1_);
    std::any b_spline_2_any = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_2_);
    std::any b_spline_3_any = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_3_);
    std::any nurbs_any = std::make_any<std::shared_ptr<spl::NURBS<1>>>(nurbs_);
    splines_ = {b_spline_1_any, b_spline_2_any, b_spline_3_any, nurbs_any};
  }

 protected:
  std::unique_ptr<io::XMLWriter> xml_writer_;
  std::vector<std::any> splines_;
};

TEST_F(A1DXMLWriter, IsCreated) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "1d_splines.xml");
  std::ifstream newFile;
  newFile.open("1d_splines.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("1d_splines.xml");
}

TEST_F(A1DXMLWriter, CreatesCorrectXMLFile) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "1d_splines.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("1d_splines.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("1d_splines.xml");
}

TEST_F(A1DXMLWriter, CreatesSplineListWith4Entries) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "1d_splines.xml");
  pugi::xml_document doc;
  doc.load_file("1d_splines.xml");
  ASSERT_STREQ(doc.first_child().name(), "SplineList");
  ASSERT_STREQ(doc.first_child().attribute("NumberOfSplines").value(), "4");
  remove("1d_splines.xml");
}

TEST_F(A1DXMLWriter, Creates4SplineEntries) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "1d_splines.xml");
  pugi::xml_document doc;
  doc.load_file("1d_splines.xml");
  pugi::xml_node spline_node = doc.child("SplineList").first_child();
  for (int i = 0; i < 4; i++, spline_node = spline_node.next_sibling()) {
    ASSERT_STREQ(spline_node.name(), "SplineEntry");
  }
  remove("1d_splines.xml");
}

TEST_F(A1DXMLWriter, CreatesNoWeightsForBSplines) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "1d_splines.xml");
  pugi::xml_document doc;
  doc.load_file("1d_splines.xml");
  pugi::xml_node spline_node = doc.child("SplineList").child("SplineEntry");
  for (int i = 0; i < 4; i++, spline_node = spline_node.next_sibling()) {
    if (i < 3) {
      ASSERT_STREQ(spline_node.child("wght").name(), "");
    } else {
      ASSERT_STREQ(spline_node.child("wght").name(), "wght");
    }
  }
  remove("1d_splines.xml");
}

TEST_F(A1DXMLWriter, WritesSplineDimension1ForAllSplines) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "1d_splines.xml");
  pugi::xml_document doc;
  doc.load_file("1d_splines.xml");
  pugi::xml_node spline_node = doc.child("SplineList").child("SplineEntry");
  for (int i = 0; i < 4; i++, spline_node = spline_node.next_sibling()) {
    ASSERT_THAT(strtod(spline_node.attribute("splDim").value(), nullptr), 1);
  }
  remove("1d_splines.xml");
}

TEST_F(A1DXMLWriter, WritesCorrectSpaceDimensions) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "1d_splines.xml");
  pugi::xml_document doc;
  doc.load_file("1d_splines.xml");
  std::vector<int> space_dimensions = {1, 2, 3, 2};
  pugi::xml_node spline_node = doc.child("SplineList").child("SplineEntry");
  for (int i = 0; i < 4; i++, spline_node = spline_node.next_sibling()) {
    ASSERT_THAT(strtod(spline_node.attribute("spaceDim").value(), nullptr), space_dimensions[i]);
  }
  remove("1d_splines.xml");
}

TEST_F(A1DXMLWriter, ReturnsSameValuesBeforeAndAfterWritingAndReadingXMLFile) {  // NOLINT
  xml_writer_->WriteXMLFile(splines_, "1d_splines.xml");
  std::unique_ptr<io::XMLReader> xml_reader(std::make_unique<io::XMLReader>());
  auto bspline_1_after = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(xml_reader->ReadXMLFile("1d_splines.xml")[0]);
  auto bspline_2_after = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(xml_reader->ReadXMLFile("1d_splines.xml")[1]);
  auto bspline_3_after = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(xml_reader->ReadXMLFile("1d_splines.xml")[2]);
  auto nurbs_after = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(xml_reader->ReadXMLFile("1d_splines.xml")[3]);

  ASSERT_THAT(b_spline_1_->Evaluate({ParamCoord(0.75839)}, {0})[0],
              DoubleEq(bspline_1_after->Evaluate({ParamCoord(0.75839)}, {0})[0]));

  ASSERT_THAT(b_spline_2_->Evaluate({ParamCoord(0.47681)}, {0})[0],
              DoubleEq(bspline_2_after->Evaluate({ParamCoord(0.47681)}, {0})[0]));
  ASSERT_THAT(b_spline_2_->Evaluate({ParamCoord(0.47681)}, {1})[0],
              DoubleEq(bspline_2_after->Evaluate({ParamCoord(0.47681)}, {1})[0]));

  ASSERT_THAT(b_spline_3_->Evaluate({ParamCoord(0.89463)}, {0})[0],
              DoubleEq(bspline_3_after->Evaluate({ParamCoord(0.89463)}, {0})[0]));
  ASSERT_THAT(b_spline_3_->Evaluate({ParamCoord(0.89463)}, {1})[0],
              DoubleEq(bspline_3_after->Evaluate({ParamCoord(0.89463)}, {1})[0]));
  ASSERT_THAT(b_spline_3_->Evaluate({ParamCoord(0.89463)}, {2})[0],
              DoubleEq(bspline_3_after->Evaluate({ParamCoord(0.89463)}, {2})[0]));

  ASSERT_THAT(nurbs_->Evaluate({ParamCoord(0.13697)}, {0})[0],
              DoubleEq(nurbs_after->Evaluate({ParamCoord(0.13697)}, {0})[0]));
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord(0.13697)}, {1})[0],
              DoubleEq(nurbs_after->Evaluate({ParamCoord(0.13697)}, {1})[0]));
  remove("1d_splines.xml");
}
