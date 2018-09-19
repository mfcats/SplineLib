/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <config.h>

#include "gmock/gmock.h"

#include "b_spline.h"
#include "irit_reader.h"

using testing::Test;
using testing::DoubleEq;

class A1DIRITReader : public Test {
 public:
  A1DIRITReader() : irit_reader(std::make_unique<io::IRITReader<1>>()) {}

 protected:
  std::unique_ptr<io::IRITReader<1>> irit_reader;
};

TEST_F(A1DIRITReader, Finds3SplinesOfDimension1) {  // NOLINT
  ASSERT_THAT(irit_reader->ReadIRITFile(path_to_iris_file).size(), 3);
}

class A1DBSplineFromIRITFile : public A1DIRITReader {
 public:
  A1DBSplineFromIRITFile() {
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

    b_spline_ = std::make_unique<spl::BSpline<1>>(knot_vector, degree, control_points);
  }

 protected:
  std::unique_ptr<spl::BSpline<1>> b_spline_;
};

TEST_F(A1DBSplineFromIRITFile, ReturnsDegree3) {  // NOLINT
  ASSERT_THAT(std::any_cast<spl::BSpline<1>>(irit_reader->ReadIRITFile(path_to_iris_file)[0]).GetDegree(0).get(), 3);
}

TEST_F(A1DBSplineFromIRITFile, ReturnsSameValueAsSplineFromIRITFile) {  // NOLINT
  std::any spline_from_file = irit_reader->ReadIRITFile(path_to_iris_file)[0];
  ASSERT_THAT(b_spline_->Evaluate({ParamCoord{0.5}}, {0})[0],
              DoubleEq(std::any_cast<spl::BSpline<1>>(spline_from_file).Evaluate({ParamCoord{0.5}}, {0})[0]));
  ASSERT_THAT(b_spline_->Evaluate({ParamCoord{10.5}}, {0})[0],
              DoubleEq(std::any_cast<spl::BSpline<1>>(spline_from_file).Evaluate({ParamCoord{10.5}}, {0})[0]));
}

class ASecond1DBSplineFromIRITFile : public A1DIRITReader {
 public:
  ASecond1DBSplineFromIRITFile() {
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

    b_spline_ = std::make_unique<spl::BSpline<1>>(knot_vector, degree, control_points);
  }
 protected:
  std::unique_ptr<spl::BSpline<1>> b_spline_;
};

TEST_F(ASecond1DBSplineFromIRITFile, ReturnsDegree2) {  // NOLINT
  ASSERT_THAT(std::any_cast<spl::BSpline<1>>(irit_reader->ReadIRITFile(path_to_iris_file)[1]).GetDegree(0).get(), 2);
}

TEST_F(ASecond1DBSplineFromIRITFile, ReturnsSameValueAsSplineFromIRITFile) {  // NOLINT
  std::any spline_from_file = irit_reader->ReadIRITFile(path_to_iris_file)[1];
  ASSERT_THAT(b_spline_->Evaluate({ParamCoord{0.5}}, {0})[0],
              DoubleEq(std::any_cast<spl::BSpline<1>>(spline_from_file).Evaluate({ParamCoord{0.5}}, {0})[0]));
  ASSERT_THAT(b_spline_->Evaluate({ParamCoord{1}}, {0})[0],
              DoubleEq(std::any_cast<spl::BSpline<1>>(spline_from_file).Evaluate({ParamCoord{1}}, {0})[0]));
}

class A1DNURBSFromIRITFile : public A1DIRITReader {
 public:
  A1DNURBSFromIRITFile() {
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

    nurbs_ = std::make_unique<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
  }
 protected:
  std::unique_ptr<spl::NURBS<1>> nurbs_;
};

TEST_F(A1DNURBSFromIRITFile, ReturnsDegree2) {  // NOLINT
  ASSERT_THAT(std::any_cast<spl::NURBS<1>>(irit_reader->ReadIRITFile(path_to_iris_file)[2]).GetDegree(0).get(), 2);
}

TEST_F(A1DNURBSFromIRITFile, ReturnsSameValueAsSplineFromIRITFile) {  // NOLINT
  std::any spline_from_file = irit_reader->ReadIRITFile(path_to_iris_file)[2];
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{0.5}}, {0})[0],
              DoubleEq(std::any_cast<spl::NURBS<1>>(spline_from_file).Evaluate({ParamCoord{0.5}}, {0})[0]));
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{1}}, {0})[0],
              DoubleEq(std::any_cast<spl::NURBS<1>>(spline_from_file).Evaluate({ParamCoord{1}}, {0})[0]));
}

class A2DIRITReader : public Test {
 public:
  A2DIRITReader() : irit_reader(std::make_unique<io::IRITReader<2>>()) {}

 protected:
  std::unique_ptr<io::IRITReader<2>> irit_reader;
};

TEST_F(A2DIRITReader, Finds2SplinesOfDimension2) {  // NOLINT
  ASSERT_THAT(irit_reader->ReadIRITFile(path_to_iris_file).size(), 2);
}

class A2DBSplineFromIRITFile : public A2DIRITReader {
 public:
  A2DBSplineFromIRITFile() {
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
    b_spline_ = std::make_unique<spl::BSpline<2>>(knot_vector, degree, control_points);
  }

 protected:
  std::unique_ptr<spl::BSpline<2>> b_spline_;
};

TEST_F(A2DBSplineFromIRITFile, ReturnsDegree2) {  // NOLINT
  ASSERT_THAT(std::any_cast<spl::BSpline<2>>(irit_reader->ReadIRITFile(path_to_iris_file)[0]).GetDegree(0).get(), 2);
  ASSERT_THAT(std::any_cast<spl::BSpline<2>>(irit_reader->ReadIRITFile(path_to_iris_file)[0]).GetDegree(1).get(), 2);
}

TEST_F(A2DBSplineFromIRITFile, ReturnsSameValueAsSplineFromIRITFile) {  // NOLINT
  std::any spline_from_file = irit_reader->ReadIRITFile(path_to_iris_file)[0];
  ASSERT_THAT(b_spline_->Evaluate({ParamCoord{0}, ParamCoord{0}}, {0})[0], DoubleEq(
      std::any_cast<spl::BSpline<2>>(spline_from_file).Evaluate({ParamCoord{0}, ParamCoord{0}}, {0})[0]));
  ASSERT_THAT(b_spline_->Evaluate({ParamCoord{0}, ParamCoord{0}}, {1})[0], DoubleEq(
      std::any_cast<spl::BSpline<2>>(spline_from_file).Evaluate({ParamCoord{0}, ParamCoord{0}}, {1})[0]));
  ASSERT_THAT(b_spline_->Evaluate({ParamCoord{0.3}, ParamCoord{0.9}}, {0})[0], DoubleEq(
      std::any_cast<spl::BSpline<2>>(spline_from_file).Evaluate({ParamCoord{0.3}, ParamCoord{0.9}}, {0})[0]));
  ASSERT_THAT(b_spline_->Evaluate({ParamCoord{0.3}, ParamCoord{0.9}}, {1})[0], DoubleEq(
      std::any_cast<spl::BSpline<2>>(spline_from_file).Evaluate({ParamCoord{0.3}, ParamCoord{0.9}}, {1})[0]));
}

class A2DNURBSFromIRITFile : public A2DIRITReader {
 public:
  A2DNURBSFromIRITFile() {
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
    nurbs_ = std::make_unique<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
  }
 protected:
  std::unique_ptr<spl::NURBS<2>> nurbs_;
};

TEST_F(A2DNURBSFromIRITFile, ReturnsDegree2) {  // NOLINT
  ASSERT_THAT(std::any_cast<spl::NURBS<2>>(irit_reader->ReadIRITFile(path_to_iris_file)[1]).GetDegree(0).get(), 2);
  ASSERT_THAT(std::any_cast<spl::NURBS<2>>(irit_reader->ReadIRITFile(path_to_iris_file)[1]).GetDegree(1).get(), 2);
}

TEST_F(A2DNURBSFromIRITFile, ReturnsSameValueAsSplineFromIRITFile) {  // NOLINT
  std::any spline_from_file = irit_reader->ReadIRITFile(path_to_iris_file)[1];
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{0}, ParamCoord{0}}, {0})[0], DoubleEq(
      std::any_cast<spl::NURBS<2>>(spline_from_file).Evaluate({ParamCoord{0}, ParamCoord{0}}, {0})[0]));
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{0}, ParamCoord{0}}, {1})[0], DoubleEq(
      std::any_cast<spl::NURBS<2>>(spline_from_file).Evaluate({ParamCoord{0}, ParamCoord{0}}, {1})[0]));
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{0.3}, ParamCoord{0.9}}, {0})[0], DoubleEq(
      std::any_cast<spl::NURBS<2>>(spline_from_file).Evaluate({ParamCoord{0.3}, ParamCoord{0.9}}, {0})[0]));
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord{0.3}, ParamCoord{0.9}}, {1})[0], DoubleEq(
      std::any_cast<spl::NURBS<2>>(spline_from_file).Evaluate({ParamCoord{0.3}, ParamCoord{0.9}}, {1})[0]));
}
