/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "src/io/xml_writer.h"

#include <fstream>

#include "gmock/gmock.h"

#include "src/io/xml_reader.h"

using testing::Test;
using testing::DoubleEq;

using namespace splinelib::src;

class A1DBSplineForXML {  // NOLINT
 public:
  A1DBSplineForXML() {
    std::array<Degree, 1> degree = {Degree{1}};
    baf::KnotVectors<1> knot_vector = {std::make_shared<baf::KnotVector>(baf::KnotVector(
        {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.5}, ParametricCoordinate{1},
         ParametricCoordinate{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0})),
        baf::ControlPoint(std::vector<double>({0.5})),
        baf::ControlPoint(std::vector<double>({1.0}))
    };
    b_spline_1d_ = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> b_spline_1d_;
  virtual ~A1DBSplineForXML() = default;
};

class A1DNURBSForXML {  // NOLINT
 public:
  A1DNURBSForXML() {
    std::array<Degree, 1> degree = {Degree{1}};
    baf::KnotVectors<1> knot_vector = {std::make_shared<baf::KnotVector>(baf::KnotVector(
        {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.5}, ParametricCoordinate{1},
         ParametricCoordinate{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0}))
    };
    std::vector<double> weights = {2.0, 1.75, 0.36};
    nurbs_1d_ = std::make_shared<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_;
  virtual ~A1DNURBSForXML() = default;
};

class A2DBSplineForXML {  // NOLINT
 public:
  A2DBSplineForXML() {
    std::array<Degree, 2> degree = {Degree{1}, Degree{2}};
    baf::KnotVectors<2> knot_vector = {
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.5}, ParametricCoordinate{1},
             ParametricCoordinate{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1},
             ParametricCoordinate{1},
             ParametricCoordinate{1}}))};
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
    b_spline_2d_ = std::make_shared<spl::BSpline<2>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> b_spline_2d_;
  virtual ~A2DBSplineForXML() = default;
};

class A2DNURBSForXML {  // NOLINT
 public:
  A2DNURBSForXML() {
    std::array<Degree, 2> degree = {Degree{1}, Degree{2}};
    baf::KnotVectors<2> knot_vector = {
        std::make_unique<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.5}, ParametricCoordinate{1},
             ParametricCoordinate{1}})),
        std::make_unique<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1},
             ParametricCoordinate{1},
             ParametricCoordinate{1}}))};
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

    nurbs_2d_ = std::make_shared<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<2>> nurbs_2d_;
  virtual ~A2DNURBSForXML() = default;
};

class A3DBSplineForXML {  // NOLINT
 public:
  A3DBSplineForXML() {
    std::array<Degree, 3> degree = {Degree{1}, Degree{1}, Degree{2}};
    baf::KnotVectors<3> knot_vector = {
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1}, ParametricCoordinate{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1}, ParametricCoordinate{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1},
             ParametricCoordinate{1},
             ParametricCoordinate{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.1, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.2, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.4, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.6, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.7, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.75, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.8, 2.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.9, 2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({0.99, 2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({0.0, 3.0, 3.0})),
        baf::ControlPoint(std::vector<double>({1.0, 3.0, 4.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0, 4.0}))
    };
    b_spline_3d_ = std::make_shared<spl::BSpline<3>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> b_spline_3d_;
  virtual ~A3DBSplineForXML() = default;
};

class A3DNURBSForXML {  // NOLINT
 public:
  A3DNURBSForXML() {
    std::array<Degree, 3> degree = {Degree{1}, Degree{2}, Degree{1}};
    baf::KnotVectors<3> knot_vector =
        {std::make_unique<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1}, ParametricCoordinate{1}})),
         std::make_unique<baf::KnotVector>(baf::KnotVector(
             {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1},
              ParametricCoordinate{1},
              ParametricCoordinate{1}})),
         std::make_unique<baf::KnotVector>(baf::KnotVector(
             {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1}, ParametricCoordinate{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.3, 0.0, 0.3, 0.0})),
        baf::ControlPoint(std::vector<double>({0.6, 0.0, 0.6, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.3, 1.0, 1.3, 1.0})),
        baf::ControlPoint(std::vector<double>({1.6, 1.0, 1.6, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0, 2.0, 1.0})),
        baf::ControlPoint(std::vector<double>({3.3, 0.6, 3.3, 0.6})),
        baf::ControlPoint(std::vector<double>({3.5, 0.8, 3.5, 0.8})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0, 4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({4.2, 0.6, 4.2, 0.6})),
        baf::ControlPoint(std::vector<double>({4.7, 0.8, 4.7, 0.8})),
        baf::ControlPoint(std::vector<double>({5.0, 1.0, 5.0, 1.0}))
    };
    std::vector<double> weights = {2.0, 1.75, 0.36, 1.0, 1.0, 0.05, 1.0, 1.0, 1.0, 2.3, 1.9, 1.4};

    nurbs_3d_ = std::make_shared<spl::NURBS<3>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<3>> nurbs_3d_;
  virtual ~A3DNURBSForXML() = default;
};

class A4DNURBSForXML {  // NOLINT
 public:
  A4DNURBSForXML() {
    std::array<Degree, 4> degree = {Degree{1}, Degree{1}, Degree{1}, Degree{1}};
    baf::KnotVectors<4> knot_vector =
        {std::make_unique<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1}, ParametricCoordinate{1}})),
         std::make_unique<baf::KnotVector>(baf::KnotVector(
             {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1}, ParametricCoordinate{1}})),
         std::make_unique<baf::KnotVector>(baf::KnotVector(
             {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1}, ParametricCoordinate{1}})),
         std::make_unique<baf::KnotVector>(baf::KnotVector(
             {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1}, ParametricCoordinate{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.1, 0.0, 0.3, 0.0})),
        baf::ControlPoint(std::vector<double>({0.2, 0.0, 0.6, 0.0})),
        baf::ControlPoint(std::vector<double>({0.5, 0.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.3, 1.0, 1.3, 1.0})),
        baf::ControlPoint(std::vector<double>({0.4, 1.0, 1.6, 1.0})),
        baf::ControlPoint(std::vector<double>({0.1, 1.0, 2.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.6, 3.3, 0.6})),
        baf::ControlPoint(std::vector<double>({-0.1, 0.8, 3.5, 0.8})),
        baf::ControlPoint(std::vector<double>({-0.3, 1.0, 4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.6, 4.2, 0.6})),
        baf::ControlPoint(std::vector<double>({0.3, 0.8, 4.7, 0.8})),
        baf::ControlPoint(std::vector<double>({0.5, 1.0, 5.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.8, 1.0, 4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.6, 4.2, 0.6})),
        baf::ControlPoint(std::vector<double>({0.4, 0.8, 4.7, 0.8})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 5.0, 1.0}))
    };
    std::vector<double> weights = {2.0, 1.75, 0.36, 1.0, 1.0, 0.05, 1.0, 1.0, 1.0, 2.3, 1.9, 1.4, 8.0, 0.01, 1.0, 3.2};

    nurbs_4d_ = std::make_shared<spl::NURBS<4>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<4>> nurbs_4d_;
  virtual ~A4DNURBSForXML() = default;
};

class AnXMLWriter : public Test, public A1DBSplineForXML, public A1DNURBSForXML, public A2DBSplineForXML,
                    public A2DNURBSForXML, public A3DBSplineForXML, public A3DNURBSForXML, public A4DNURBSForXML {
 public:
  AnXMLWriter() : xml_writer_(std::make_unique<io::XMLWriter>()) {
    std::any b_spline_1d_any = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_1d_);
    std::any nurbs_1d_any = std::make_any<std::shared_ptr<spl::NURBS<1>>>(nurbs_1d_);
    std::any b_spline_2d_any = std::make_any<std::shared_ptr<spl::BSpline<2>>>(b_spline_2d_);
    std::any nurbs_2d_any = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbs_2d_);
    std::any b_spline_3d_any = std::make_any<std::shared_ptr<spl::BSpline<3>>>(b_spline_3d_);
    std::any nurbs_3d_any = std::make_any<std::shared_ptr<spl::NURBS<3>>>(nurbs_3d_);
    std::any nurbs_4d_any = std::make_any<std::shared_ptr<spl::NURBS<4>>>(nurbs_4d_);
    splines_ =
        {b_spline_1d_any, nurbs_1d_any, b_spline_2d_any, nurbs_2d_any, b_spline_3d_any, nurbs_3d_any, nurbs_4d_any};
  }

 protected:
  std::unique_ptr<io::XMLWriter> xml_writer_;
  std::vector<std::any> splines_;
};

TEST_F(AnXMLWriter, IsCreated) {  // NOLINT
  xml_writer_->WriteFile(splines_, "splines.xml");
  std::ifstream newFile;
  newFile.open("splines.xml");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("splines.xml");
}

TEST_F(AnXMLWriter, CreatesCorrectXMLFile) {  // NOLINT
  xml_writer_->WriteFile(splines_, "splines.xml");
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("splines.xml");
  ASSERT_STREQ(result.description(), "No error");
  remove("splines.xml");
}

TEST_F(AnXMLWriter, CreatesSplineListWith7Entries) {  // NOLINT
  xml_writer_->WriteFile(splines_, "splines.xml");
  pugi::xml_document doc;
  doc.load_file("splines.xml");
  ASSERT_STREQ(doc.first_child().name(), "SplineList");
  ASSERT_STREQ(doc.first_child().attribute("NumberOfSplines").value(), "7");
  remove("splines.xml");
}

TEST_F(AnXMLWriter, Creates7SplineEntries) {  // NOLINT
  xml_writer_->WriteFile(splines_, "splines.xml");
  pugi::xml_document doc;
  doc.load_file("splines.xml");
  pugi::xml_node spline_node = doc.child("SplineList").first_child();
  for (int i = 0; i < 7; i++, spline_node = spline_node.next_sibling()) {
    ASSERT_STREQ(spline_node.name(), "SplineEntry");
  }
  remove("splines.xml");
}

TEST_F(AnXMLWriter, CreatesNoWeightsForBSplines) {  // NOLINT
  xml_writer_->WriteFile(splines_, "splines.xml");
  pugi::xml_document doc;
  doc.load_file("splines.xml");
  pugi::xml_node spline_node = doc.child("SplineList").child("SplineEntry");
  for (int i = 0; i < 7; i++, spline_node = spline_node.next_sibling()) {
    if (i % 2 == 0 && i != 6) {
      ASSERT_STREQ(spline_node.child("wght").name(), "");
    } else {
      ASSERT_STREQ(spline_node.child("wght").name(), "wght");
    }
  }
  remove("splines.xml");
}

TEST_F(AnXMLWriter, WritesCorrestSplineDimensions) {  // NOLINT
  xml_writer_->WriteFile(splines_, "splines.xml");
  pugi::xml_document doc;
  doc.load_file("splines.xml");
  pugi::xml_node spline_node = doc.child("SplineList").child("SplineEntry");
  for (int i = 0; i < 7; i++, spline_node = spline_node.next_sibling()) {
    ASSERT_THAT(strtod(spline_node.attribute("splDim").value(), nullptr), ceil((i + 1) / 2.0));
  }
  remove("splines.xml");
}

TEST_F(AnXMLWriter, WritesCorrectSpaceDimensions) {  // NOLINT
  xml_writer_->WriteFile(splines_, "splines.xml");
  pugi::xml_document doc;
  doc.load_file("splines.xml");
  std::vector<int> space_dimensions = {1, 2, 3, 2, 3, 4, 4};
  pugi::xml_node spline_node = doc.child("SplineList").child("SplineEntry");
  for (int i = 0; i < 7; i++, spline_node = spline_node.next_sibling()) {
    ASSERT_THAT(strtod(spline_node.attribute("spaceDim").value(), nullptr), space_dimensions[i]);
  }
  remove("splines.xml");
}

TEST_F(AnXMLWriter, ReturnsSameValuesBeforeAndAfterWritingAndReadingXMLFile) {  // NOLINT
  xml_writer_->WriteFile(splines_, "splines.xml");
  std::unique_ptr<io::XMLReader> xml_reader(std::make_unique<io::XMLReader>());
  auto bspline_1d_after = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(xml_reader->ReadFile("splines.xml")[0]);
  auto nurbs_1d_after = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(xml_reader->ReadFile("splines.xml")[1]);
  auto bspline_2d_after = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(xml_reader->ReadFile("splines.xml")[2]);
  auto nurbs_2d_after = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(xml_reader->ReadFile("splines.xml")[3]);
  auto bspline_3d_after = std::any_cast<std::shared_ptr<spl::BSpline<3>>>(xml_reader->ReadFile("splines.xml")[4]);
  auto nurbs_3d_after = std::any_cast<std::shared_ptr<spl::NURBS<3>>>(xml_reader->ReadFile("splines.xml")[5]);
  auto nurbs_4d_after = std::any_cast<std::shared_ptr<spl::NURBS<4>>>(xml_reader->ReadFile("splines.xml")[6]);
  ASSERT_THAT(b_spline_1d_->Evaluate({ParametricCoordinate(0.75839)}, {0})[0],
              DoubleEq(bspline_1d_after->Evaluate({ParametricCoordinate(0.75839)}, {0})[0]));

  ASSERT_THAT(nurbs_1d_->Evaluate({ParametricCoordinate(0.13697)}, {0})[0],
              DoubleEq(nurbs_1d_after->Evaluate({ParametricCoordinate(0.13697)}, {0})[0]));
  ASSERT_THAT(nurbs_1d_->Evaluate({ParametricCoordinate(0.13697)}, {1})[0],
              DoubleEq(nurbs_1d_after->Evaluate({ParametricCoordinate(0.13697)}, {1})[0]));

  ASSERT_THAT(b_spline_2d_->Evaluate({ParametricCoordinate(0.47681)}, {0})[0],
              DoubleEq(bspline_2d_after->Evaluate({ParametricCoordinate(0.47681)}, {0})[0]));
  ASSERT_THAT(b_spline_2d_->Evaluate({ParametricCoordinate(0.47681)}, {1})[0],
              DoubleEq(bspline_2d_after->Evaluate({ParametricCoordinate(0.47681)}, {1})[0]));
  ASSERT_THAT(b_spline_2d_->Evaluate({ParametricCoordinate(0.47681)}, {2})[0],
              DoubleEq(bspline_2d_after->Evaluate({ParametricCoordinate(0.47681)}, {2})[0]));

  ASSERT_THAT(nurbs_2d_->Evaluate({ParametricCoordinate(0.27856)}, {0})[0],
              DoubleEq(nurbs_2d_after->Evaluate({ParametricCoordinate(0.27856)}, {0})[0]));
  ASSERT_THAT(nurbs_2d_->Evaluate({ParametricCoordinate(0.27856)}, {1})[0],
              DoubleEq(nurbs_2d_after->Evaluate({ParametricCoordinate(0.27856)}, {1})[0]));

  ASSERT_THAT(b_spline_3d_->Evaluate({ParametricCoordinate(0.78781)}, {0})[0],
              DoubleEq(bspline_3d_after->Evaluate({ParametricCoordinate(0.78781)}, {0})[0]));
  ASSERT_THAT(b_spline_3d_->Evaluate({ParametricCoordinate(0.78781)}, {1})[0],
              DoubleEq(bspline_3d_after->Evaluate({ParametricCoordinate(0.78781)}, {1})[0]));
  ASSERT_THAT(b_spline_3d_->Evaluate({ParametricCoordinate(0.78781)}, {2})[0],
              DoubleEq(bspline_3d_after->Evaluate({ParametricCoordinate(0.78781)}, {2})[0]));

  ASSERT_THAT(nurbs_3d_->Evaluate({ParametricCoordinate(0.90069)}, {0})[0],
              DoubleEq(nurbs_3d_after->Evaluate({ParametricCoordinate(0.90069)}, {0})[0]));
  ASSERT_THAT(nurbs_3d_->Evaluate({ParametricCoordinate(0.90069)}, {1})[0],
              DoubleEq(nurbs_3d_after->Evaluate({ParametricCoordinate(0.90069)}, {1})[0]));
  ASSERT_THAT(nurbs_3d_->Evaluate({ParametricCoordinate(0.90069)}, {2})[0],
              DoubleEq(nurbs_3d_after->Evaluate({ParametricCoordinate(0.90069)}, {2})[0]));
  ASSERT_THAT(nurbs_3d_->Evaluate({ParametricCoordinate(0.90069)}, {3})[0],
              DoubleEq(nurbs_3d_after->Evaluate({ParametricCoordinate(0.90069)}, {3})[0]));

  ASSERT_THAT(nurbs_4d_->Evaluate({ParametricCoordinate(0.3574)}, {0})[0],
              DoubleEq(nurbs_4d_after->Evaluate({ParametricCoordinate(0.3574)}, {0})[0]));
  ASSERT_THAT(nurbs_4d_->Evaluate({ParametricCoordinate(0.3574)}, {1})[0],
              DoubleEq(nurbs_4d_after->Evaluate({ParametricCoordinate(0.3574)}, {1})[0]));
  ASSERT_THAT(nurbs_4d_->Evaluate({ParametricCoordinate(0.3574)}, {2})[0],
              DoubleEq(nurbs_4d_after->Evaluate({ParametricCoordinate(0.3574)}, {2})[0]));
  ASSERT_THAT(nurbs_4d_->Evaluate({ParametricCoordinate(0.3574)}, {3})[0],
              DoubleEq(nurbs_4d_after->Evaluate({ParametricCoordinate(0.3574)}, {3})[0]));
  remove("splines.xml");
}

TEST_F(AnXMLWriter, ThrowsForSplineOfDimensionFive) {  // NOLINT
  std::shared_ptr<spl::NURBS<5>> nurbs_5d_;
  std::any nurbs_5d_any = std::make_any<std::shared_ptr<spl::NURBS<5>>>(nurbs_5d_);
  ASSERT_THROW(xml_writer_->WriteFile({nurbs_5d_any}, "5d_spline.xml"), std::runtime_error);
}
