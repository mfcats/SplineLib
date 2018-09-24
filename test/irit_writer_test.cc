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

#include "b_spline.h"
#include "irit_writer.h"
#include "nurbs.h"

using testing::Test;
using testing::Ne;

class A1DIRITWriter : public Test {
 public:
  A1DIRITWriter() {
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1},
                         ParamCoord{1}, ParamCoord{2}, ParamCoord{3}, ParamCoord{4}, ParamCoord{5}, ParamCoord{6},
                         ParamCoord{7}, ParamCoord{8}, ParamCoord{9}, ParamCoord{10}, ParamCoord{11}, ParamCoord{11},
                         ParamCoord{11}, ParamCoord{11}}))};
    std::array<Degree, 1> degree = {Degree{3}};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.874, 0})),
        baf::ControlPoint(std::vector<double>({0.899333, 0.0253333})),
        baf::ControlPoint(std::vector<double>({0.924667, 0.0506667})),
        baf::ControlPoint(std::vector<double>({0.95, 0.076})),
        baf::ControlPoint(std::vector<double>({0.95, 0.76})),
        baf::ControlPoint(std::vector<double>({0.304, 1.52})),
        baf::ControlPoint(std::vector<double>({0.304, 1.9})),
        baf::ControlPoint(std::vector<double>({0.494, 2.09})),
        baf::ControlPoint(std::vector<double>({0.722, 2.242})),
        baf::ControlPoint(std::vector<double>({0.722, 2.318})),
        baf::ControlPoint(std::vector<double>({0.38, 2.508})),
        baf::ControlPoint(std::vector<double>({0.418, 2.698})),
        baf::ControlPoint(std::vector<double>({0.57, 2.812})),
        baf::ControlPoint(std::vector<double>({0.57, 3.42})),
        baf::ControlPoint(std::vector<double>({0.19, 3.572})),
        baf::ControlPoint(std::vector<double>({0, 3.572}))
    };
    std::shared_ptr<spl::BSpline<1>>
        b_spline_ptr = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
    std::any b_spline_ = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_ptr);
    splines.push_back(b_spline_);
    irit_writer = std::make_unique<io::IRITWriter<1>>(splines);
  }

 protected:
  std::unique_ptr<io::IRITWriter<1>> irit_writer;
  std::vector<std::any> splines;
};

TEST_F(A1DIRITWriter, IsCreated) {  // NOLINT
  irit_writer->WriteIRITFile("bspline.itd");
  std::ifstream newFile;
  newFile.open("bspline.itd");
  ASSERT_TRUE(newFile.is_open());
  newFile.close();
  remove("bspline.itd");
}

TEST_F(A1DIRITWriter, CreatesCorrectFile) {  // NOLINT
  irit_writer->WriteIRITFile("bspline.itd");
  std::ifstream newFile;
  newFile.open("bspline.itd");
  std::string line, file;
  while (getline(newFile, line)) {
    file += line;
  }
  ASSERT_THAT(file.find("BSPLINE 16 4"), Ne(std::string::npos));
  ASSERT_THAT(file.find("CURVE"), Ne(std::string::npos));
  ASSERT_THAT(file.find("E2"), Ne(std::string::npos));
  ASSERT_THAT(file.find("KV"), Ne(std::string::npos));
  ASSERT_THAT(file.find("11.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("3.572000]"), Ne(std::string::npos));
  // remove("bspline.itd");
}

class A1DIRITWriterWithTwoSplines : public A1DIRITWriter {
 public:
  A1DIRITWriterWithTwoSplines() {
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.16666666666667},
                         ParamCoord{0.33333333333333}, ParamCoord{0.5}, ParamCoord{0.66666666666667},
                         ParamCoord{0.83333333333333}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    std::array<Degree, 1> degree = {Degree{2}};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0, 0, 0})),
        baf::ControlPoint(std::vector<double>({0.5, 0.5, 0})),
        baf::ControlPoint(std::vector<double>({0.5, 0.5, 1})),
        baf::ControlPoint(std::vector<double>({1, 1, 1})),
        baf::ControlPoint(std::vector<double>({2, 1, 0})),
        baf::ControlPoint(std::vector<double>({1.5, 0.5, -0.5})),
        baf::ControlPoint(std::vector<double>({0.8, 0.3, -0.4})),
        baf::ControlPoint(std::vector<double>({0.5, 0, 0}))
    };
    std::shared_ptr<spl::BSpline<1>>
        b_spline_ptr = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
    std::any b_spline_ = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_ptr);
    splines.push_back(b_spline_);
    irit_writer = std::make_unique<io::IRITWriter<1>>(splines);
  }
};

TEST_F(A1DIRITWriterWithTwoSplines, CreatesCorrectFile) {  // NOLINT
  irit_writer->WriteIRITFile("2bsplines.itd");
  std::ifstream newFile;
  newFile.open("2bsplines.itd");
  std::string line, file;
  while (getline(newFile, line)) {
    file += line;
  }
  ASSERT_THAT(file.find("BSPLINE 16 4"), Ne(std::string::npos));
  ASSERT_THAT(file.find("BSPLINE 8 3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("CURVE"), Ne(std::string::npos));
  ASSERT_THAT(file.find("E2"), Ne(std::string::npos));
  ASSERT_THAT(file.find("E3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("KV"), Ne(std::string::npos));
  ASSERT_THAT(file.find("11.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("1.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("3.572000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("-0.400000]"), Ne(std::string::npos));
  // remove("2bsplines.itd");
}

class A1DIRITWriterWithNURBS : public A1DIRITWriterWithTwoSplines {
 public:
  A1DIRITWriterWithNURBS() {
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.25}, ParamCoord{0.5},
                         ParamCoord{0.75}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    std::array<Degree, 1> degree = {Degree{2}};
    std::vector<double> weights = {0.5, 0.5, 0.5, 1, 1, 1};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0, 0})),
        baf::ControlPoint(std::vector<double>({1, 0.3})),
        baf::ControlPoint(std::vector<double>({2, 0.5})),
        baf::ControlPoint(std::vector<double>({0, 0.6})),
        baf::ControlPoint(std::vector<double>({1, 0.8})),
        baf::ControlPoint(std::vector<double>({2, 1}))
    };
    std::shared_ptr<spl::NURBS<1>>
        nurbs_ptr = std::make_shared<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
    std::any b_spline_ = std::make_any<std::shared_ptr<spl::NURBS<1>>>(nurbs_ptr);
    splines.push_back(b_spline_);
    irit_writer = std::make_unique<io::IRITWriter<1>>(splines);
  }
};

TEST_F(A1DIRITWriterWithNURBS, CreatesCorrectFile) {  // NOLINT
  irit_writer->WriteIRITFile("nurbs.itd");
  std::ifstream newFile;
  newFile.open("nurbs.itd");
  std::string line, file;
  while (getline(newFile, line)) {
    file += line;
  }
  ASSERT_THAT(file.find("BSPLINE 16 4"), Ne(std::string::npos));
  ASSERT_THAT(file.find("BSPLINE 8 3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("CURVE"), Ne(std::string::npos));
  ASSERT_THAT(file.find("E2"), Ne(std::string::npos));
  ASSERT_THAT(file.find("E3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("KV"), Ne(std::string::npos));
  ASSERT_THAT(file.find("11.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("1.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("3.572000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("-0.400000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[0.500000 2.000000 0.500000]"), Ne(std::string::npos));
  // remove("nurbs.itd");
}

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
  ASSERT_THAT(file.find("BSPLINE 3 3 3 3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("SURFACE"), Ne(std::string::npos));
  ASSERT_THAT(file.find("E3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("KV"), Ne(std::string::npos));
  ASSERT_THAT(file.find("1.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("0.500000]"), Ne(std::string::npos));
  // remove("2d_bspline.itd");
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
  ASSERT_THAT(file.find("BSPLINE 3 3 3 3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("SURFACE"), Ne(std::string::npos));
  ASSERT_THAT(file.find("E3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("KV"), Ne(std::string::npos));
  ASSERT_THAT(file.find("1.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("0.500000]"), Ne(std::string::npos));
  // remove("2d_nurbs.itd");
}
