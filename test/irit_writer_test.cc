/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "gmock/gmock.h"

#include "irit_writer.h"
#include "irit_reader.h"

using testing::Test;
using testing::Ne;
using testing::DoubleEq;


class A2DIRITWriter : public Test {
 public:
  A2DIRITWriter() {
    std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector = {
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}, ParamCoord{1}}))};
    std::array<Degree, 2> degree = {Degree{2}, Degree{2}};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0, 0, 0})),
        baf::ControlPoint(std::vector<double>({0, 1, 0.3})),
        baf::ControlPoint(std::vector<double>({0, 2, 0.5})),
        baf::ControlPoint(std::vector<double>({1, 0, 0.6})),
        baf::ControlPoint(std::vector<double>({1, 1, 0.8})),
        baf::ControlPoint(std::vector<double>({1, 2, 1})),
        baf::ControlPoint(std::vector<double>({2, 0, 0.5})),
        baf::ControlPoint(std::vector<double>({2, 1, 0.9})),
        baf::ControlPoint(std::vector<double>({2, 2, 0.5}))
    };
    std::shared_ptr<spl::BSpline<2>>
        b_spline_ptr = std::make_shared<spl::BSpline<2>>(knot_vector, degree, control_points);
    std::any b_spline_ = std::make_any<std::shared_ptr<spl::BSpline<2>>>(b_spline_ptr);
    splines.push_back(b_spline_);
    irit_writer = std::make_unique<io::IRITWriter<2>>(splines);
  }

 protected:
  std::unique_ptr<io::IRITWriter<2>> irit_writer;
  std::vector<std::any> splines;
};

TEST_F(A2DIRITWriter, CreatesCorrectFile) {  // NOLINT
  irit_writer->WriteIRITFile("2d_bspline.itd");
  std::ifstream newFile;
  newFile.open("2d_bspline.itd");
  std::string line, file;
  while (getline(newFile, line)) {
    file += line;
  }
  ASSERT_THAT(file.find("SURFACE BSPLINE 3 3 3 3 E3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[KV "), Ne(std::string::npos));
  ASSERT_THAT(file.find(" 1.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[2.000000 2.000000 0.500000]"), Ne(std::string::npos));
  remove("2d_bspline.itd");
}

class A2DIRITWriterwithNURBS : public A2DIRITWriter {
 public:
  A2DIRITWriterwithNURBS() {
    std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector = {
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}, ParamCoord{1}}))};
    std::array<Degree, 2> degree = {Degree{2}, Degree{2}};
    std::vector<double> weights = {0.5, 0.5, 0.5, 1, 1, 1, 2, 2, 2};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0, 0})),
        baf::ControlPoint(std::vector<double>({1, 0.3})),
        baf::ControlPoint(std::vector<double>({2, 0.5})),
        baf::ControlPoint(std::vector<double>({0, 0.6})),
        baf::ControlPoint(std::vector<double>({1, 0.8})),
        baf::ControlPoint(std::vector<double>({2, 1})),
        baf::ControlPoint(std::vector<double>({0, 0.5})),
        baf::ControlPoint(std::vector<double>({1, 0.9})),
        baf::ControlPoint(std::vector<double>({2, 0.5}))
    };
    std::shared_ptr<spl::NURBS<2>>
        nurbs_ptr = std::make_shared<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
    std::any nurbs = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbs_ptr);
    splines.push_back(nurbs);
    irit_writer = std::make_unique<io::IRITWriter<2>>(splines);
  }
};

TEST_F(A2DIRITWriterwithNURBS, CreatesCorrectFile) {  // NOLINT
  irit_writer->WriteIRITFile("2d_nurbs.itd");
  std::ifstream newFile;
  newFile.open("2d_nurbs.itd");
  std::string line, file;
  while (getline(newFile, line)) {
    file += line;
  }
  ASSERT_THAT(file.find("SURFACE BSPLINE 3 3 3 3 E3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("SURFACE BSPLINE 3 3 3 3 P2"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[KV "), Ne(std::string::npos));
  ASSERT_THAT(file.find(" 1.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[2.000000 2.000000 0.500000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[0.500000 0.000000 0.000000]"), Ne(std::string::npos));
  remove("2d_nurbs.itd");
}

TEST_F(A2DIRITWriterwithNURBS, ReturnsSameValuesBeforeAndAfterWritingAndReadingIRITFile) {  // NOLINT
  spl::BSpline<2> bspline_before = *std::any_cast<std::shared_ptr<spl::BSpline<2>>>(splines[0]);
  spl::NURBS<2> nurbs_before = *std::any_cast<std::shared_ptr<spl::NURBS<2>>>(splines[1]);
  irit_writer->WriteIRITFile("splines.itd");
  std::unique_ptr<io::IRITReader<2>> irit_reader(std::make_unique<io::IRITReader<2>>());
  auto bspline_after = std::any_cast<spl::BSpline<2>>(irit_reader->ReadIRITFile("splines.itd")[0]);
  auto nurbs_after = std::any_cast<spl::NURBS<2>>(irit_reader->ReadIRITFile("splines.itd")[1]);
  ASSERT_THAT(bspline_before.Evaluate({ParamCoord(0.75839), ParamCoord(0.01453)}, {0})[0],
              DoubleEq(bspline_after.Evaluate({ParamCoord(0.75839), ParamCoord(0.01453)}, {0})[0]));
  ASSERT_THAT(nurbs_before.Evaluate({ParamCoord(0.75839), ParamCoord(0.01453)}, {0})[0],
              DoubleEq(nurbs_after.Evaluate({ParamCoord(0.75839), ParamCoord(0.01453)}, {0})[0]));
  remove("splines.itd");
}

class A3DIRITWriter : public Test {
 public:
  A3DIRITWriter() {
    std::array<std::shared_ptr<baf::KnotVector>, 3> knot_vector = {
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}}))};
    std::array<Degree, 3> degree = {Degree{1}, Degree{1}, Degree{1}};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0, 0, 0})),
        baf::ControlPoint(std::vector<double>({1, 0, 0})),
        baf::ControlPoint(std::vector<double>({0, 1, 0})),
        baf::ControlPoint(std::vector<double>({1, 1, 0})),
        baf::ControlPoint(std::vector<double>({0, 0, 1})),
        baf::ControlPoint(std::vector<double>({1, 0, 1})),
        baf::ControlPoint(std::vector<double>({0, 1, 1})),
        baf::ControlPoint(std::vector<double>({1, 1, 1}))
    };
    std::shared_ptr<spl::BSpline<3>>
        b_spline_ptr = std::make_shared<spl::BSpline<3>>(knot_vector, degree, control_points);
    splines.push_back(std::make_any<std::shared_ptr<spl::BSpline<3>>>(b_spline_ptr));

    knot_vector = {
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}}))};
    degree = {Degree{1}, Degree{1}, Degree{1}};
    std::vector<double> weights = {0.2, 0.3, 0.5, 0.75, 1, 1.3, 1.5, 2};
    control_points = {
        baf::ControlPoint(std::vector<double>({0, 0, 0})),
        baf::ControlPoint(std::vector<double>({1, 0, 0})),
        baf::ControlPoint(std::vector<double>({0, 1, 0})),
        baf::ControlPoint(std::vector<double>({1, 1, 0})),
        baf::ControlPoint(std::vector<double>({0, 0, 1})),
        baf::ControlPoint(std::vector<double>({1, 0, 1})),
        baf::ControlPoint(std::vector<double>({0, 1, 1})),
        baf::ControlPoint(std::vector<double>({1, 1, 1}))
    };
    std::shared_ptr<spl::NURBS<3>>
        nurbs_ptr = std::make_shared<spl::NURBS<3>>(knot_vector, degree, control_points, weights);
    splines.push_back(std::make_any<std::shared_ptr<spl::NURBS<3>>>(nurbs_ptr));

    irit_writer = std::make_unique<io::IRITWriter<3>>(splines);
  }

 protected:
  std::unique_ptr<io::IRITWriter<3>> irit_writer;
  std::vector<std::any> splines;
};

TEST_F(A3DIRITWriter, ReturnsSameValuesBeforeAndAfterWritingAndReadingIRITFile) {  // NOLINT
  spl::BSpline<3> bspline_before = *std::any_cast<std::shared_ptr<spl::BSpline<3>>>(splines[0]);
  spl::NURBS<3> nurbs_before = *std::any_cast<std::shared_ptr<spl::NURBS<3>>>(splines[1]);
  irit_writer->WriteIRITFile("3d_splines.itd");
  std::unique_ptr<io::IRITReader<3>> irit_reader(std::make_unique<io::IRITReader<3>>());
  auto bspline_after = std::any_cast<spl::BSpline<3>>(irit_reader->ReadIRITFile("3d_splines.itd")[0]);
  auto nurbs_after = std::any_cast<spl::NURBS<3>>(irit_reader->ReadIRITFile("3d_splines.itd")[1]);
  ASSERT_THAT(bspline_before.Evaluate({ParamCoord(0.75839), ParamCoord(0.01453), ParamCoord(0.5789)}, {0})[0],
              DoubleEq(bspline_after.Evaluate({ParamCoord(0.75839), ParamCoord(0.01453), ParamCoord(0.5789)}, {0})[0]));
  ASSERT_THAT(nurbs_before.Evaluate({ParamCoord(0.75839), ParamCoord(0.01453), ParamCoord(0.5789)}, {0})[0],
              DoubleEq(nurbs_after.Evaluate({ParamCoord(0.75839), ParamCoord(0.01453), ParamCoord(0.5789)}, {0})[0]));
  remove("3d_splines.itd");
}
