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

#include "irit_reader.h"
#include "irit_writer.h"
#include "xml_reader.h"
#include "xml_writer.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;
using testing::Ne;

class A1DBSplineForIRIT {  // NOLINT
 public:
  A1DBSplineForIRIT() {
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
    b_spline_1d_ = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> b_spline_1d_;
  virtual ~A1DBSplineForIRIT() = default;
};

class A1DNURBSForIRIT {  // NOLINT
 public:
  A1DNURBSForIRIT() {
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
    nurbs_1d_ = std::make_shared<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_;
  virtual ~A1DNURBSForIRIT() = default;
};

class A2DBSplineForIRIT {  // NOLINT
 public:
  A2DBSplineForIRIT() {
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
    b_spline_2d_ = std::make_shared<spl::BSpline<2>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> b_spline_2d_;
  virtual ~A2DBSplineForIRIT() = default;
};

class A2DNURBSForIRIT {  // NOLINT
 public:
  A2DNURBSForIRIT() {
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
    nurbs_2d_ = std::make_shared<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<2>> nurbs_2d_;
  virtual ~A2DNURBSForIRIT() = default;
};

class A3DBSplineForIRIT {  // NOLINT
 public:
  A3DBSplineForIRIT() {
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
    b_spline_3d_ = std::make_shared<spl::BSpline<3>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> b_spline_3d_;
  virtual ~A3DBSplineForIRIT() = default;
};

class A3DNURBSForIRIT {  // NOLINT
 public:
  A3DNURBSForIRIT() {
    std::array<std::shared_ptr<baf::KnotVector>, 3> knot_vector = {
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}}))};
    std::array<Degree, 3> degree = {Degree{1}, Degree{1}, Degree{1}};
    std::vector<double> weights = {0.2, 0.3, 0.5, 0.75, 1, 1.3, 1.5, 2};
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
    nurbs_3d_ = std::make_shared<spl::NURBS<3>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<3>> nurbs_3d_;
  virtual ~A3DNURBSForIRIT() = default;
};

class AnIRITReader : public Test, public A1DBSplineForIRIT, public A1DNURBSForIRIT, public A2DBSplineForIRIT,
                     public A2DNURBSForIRIT, public A3DBSplineForIRIT, public A3DNURBSForIRIT {
 public:
  AnIRITReader() : irit_reader(std::make_unique<io::IRITReader>()) {}

 protected:
  std::unique_ptr<io::IRITReader> irit_reader;
};

TEST_F(AnIRITReader, ThrowsExceptionForNonExistingFile) {  // NOLINT
  ASSERT_THROW(irit_reader->ReadFile("testing.itd"), std::runtime_error);
}

TEST_F(AnIRITReader, Finds6Splines) {  // NOLINT
  ASSERT_THAT(irit_reader->ReadFile(path_to_irit_file).size(), 6);
}

TEST_F(AnIRITReader, ReturnsCorrectDegree) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<1>>>(
      irit_reader->ReadFile(path_to_irit_file)[0])->GetDegree(0).get(), b_spline_1d_->GetDegree(0).get());
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<1>>>(
      irit_reader->ReadFile(path_to_irit_file)[1])->GetDegree(0).get(), nurbs_1d_->GetDegree(0).get());
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      irit_reader->ReadFile(path_to_irit_file)[2])->GetDegree(1).get(), b_spline_2d_->GetDegree(1).get());
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      irit_reader->ReadFile(path_to_irit_file)[3])->GetDegree(1).get(), nurbs_2d_->GetDegree(1).get());
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<3>>>(
      irit_reader->ReadFile(path_to_irit_file)[4])->GetDegree(2).get(), b_spline_3d_->GetDegree(2).get());
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<3>>>(
      irit_reader->ReadFile(path_to_irit_file)[5])->GetDegree(2).get(), nurbs_3d_->GetDegree(2).get());
}

TEST_F(AnIRITReader, ReturnsSameValuesAsGivenSplines) {  // NOLINT
  std::vector<std::any> spline_vector = irit_reader->ReadFile(path_to_irit_file);
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<1>>>(spline_vector[0])->Evaluate({ParamCoord{0.5}}, {0})[0],
              DoubleEq(b_spline_1d_->Evaluate({ParamCoord{0.5}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<1>>>(spline_vector[1])->Evaluate({ParamCoord{0.123}}, {0})[0],
              DoubleEq(nurbs_1d_->Evaluate({ParamCoord{0.123}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(spline_vector[2])->Evaluate({ParamCoord{0.5}}, {1})[0],
              DoubleEq(b_spline_2d_->Evaluate({ParamCoord{0.5}}, {1})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(spline_vector[3])->Evaluate({ParamCoord{0.375}}, {1})[0],
              DoubleEq(nurbs_2d_->Evaluate({ParamCoord{0.375}}, {1})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<3>>>(spline_vector[4])->Evaluate({ParamCoord{0.999}}, {2})[0],
              DoubleEq(b_spline_3d_->Evaluate({ParamCoord{0.999}}, {2})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<3>>>(spline_vector[5])->Evaluate({ParamCoord{0.007}}, {2})[0],
              DoubleEq(nurbs_3d_->Evaluate({ParamCoord{0.007}}, {2})[0]));
}

class AnIRITWriter : public Test, public A1DBSplineForIRIT, public A1DNURBSForIRIT, public A2DBSplineForIRIT,
                     public A2DNURBSForIRIT, public A3DBSplineForIRIT, public A3DNURBSForIRIT {
 public:
  AnIRITWriter() : irit_writer_(std::make_unique<io::IRITWriter>()) {
    std::any b_spline_1d_any = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_1d_);
    std::any nurbs_1d_any = std::make_any<std::shared_ptr<spl::NURBS<1>>>(nurbs_1d_);
    std::any b_spline_2d_any = std::make_any<std::shared_ptr<spl::BSpline<2>>>(b_spline_2d_);
    std::any nurbs_2d_any = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbs_2d_);
    std::any b_spline_3d_any = std::make_any<std::shared_ptr<spl::BSpline<3>>>(b_spline_3d_);
    std::any nurbs_3d_any = std::make_any<std::shared_ptr<spl::NURBS<3>>>(nurbs_3d_);
    splines_ = {b_spline_1d_any, nurbs_1d_any, b_spline_2d_any, nurbs_2d_any, b_spline_3d_any, nurbs_3d_any};
  }

 protected:
  std::unique_ptr<io::IRITWriter> irit_writer_;
  std::vector<std::any> splines_;
};

TEST_F(AnIRITWriter, ThrowsExceptionForWrongInputType) {  // NOLINT
  ASSERT_THROW(io::IRITWriter().WriteFile({std::make_any<int>(8)}, "splines.itd"), std::runtime_error);
  remove("splines.itd");
}

TEST_F(AnIRITWriter, CreatesCorrectFile) {  // NOLINT
  irit_writer_->WriteFile(splines_, "splines.itd");
  std::ifstream newFile;
  newFile.open("splines.itd");
  std::string line, file;
  while (getline(newFile, line)) {
    file += line;
  }
  ASSERT_THAT(file.find("[CURVE BSPLINE 8 3 E3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[CURVE BSPLINE 6 3 P2"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[SURFACE BSPLINE 3 3 3 3 E3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[SURFACE BSPLINE 3 3 3 3 P2"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[TRIVAR BSPLINE 2 2 2 2 2 2 E3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[TRIVAR BSPLINE 2 2 2 2 2 2 P3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[KV "), Ne(std::string::npos));
  ASSERT_THAT(file.find(" 1.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[0.800000 0.300000 -0.400000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[0.500000 2.000000 0.500000]"), Ne(std::string::npos));
  remove("splines.itd");
}

TEST_F(AnIRITWriter, ReturnsSameValuesBeforeAndAfterWritingAndReadingIRITFile) {  // NOLINT
  irit_writer_->WriteFile(splines_, "splines.itd");
  std::unique_ptr<io::IRITReader> irit_reader(std::make_unique<io::IRITReader>());
  ASSERT_THAT(irit_reader->ReadFile("splines.itd").size(), 6);
  auto bspline_1d_after = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(irit_reader->ReadFile("splines.itd")[0]);
  auto nurbs_1d_after = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(irit_reader->ReadFile("splines.itd")[1]);
  auto bspline_2d_after = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(irit_reader->ReadFile("splines.itd")[2]);
  auto nurbs_2d_after = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_reader->ReadFile("splines.itd")[3]);
  auto bspline_3d_after = std::any_cast<std::shared_ptr<spl::BSpline<3>>>(irit_reader->ReadFile("splines.itd")[4]);
  auto nurbs_3d_after = std::any_cast<std::shared_ptr<spl::NURBS<3>>>(irit_reader->ReadFile("splines.itd")[5]);
  ASSERT_THAT(b_spline_1d_->Evaluate({ParamCoord(0.75839)}, {0})[0],
              DoubleNear(bspline_1d_after->Evaluate({ParamCoord(0.75839)}, {0})[0], 0.00001));
  ASSERT_THAT(b_spline_1d_->Evaluate({ParamCoord(0.75839)}, {1})[0],
              DoubleNear(bspline_1d_after->Evaluate({ParamCoord(0.75839)}, {1})[0], 0.00001));
  ASSERT_THAT(b_spline_1d_->Evaluate({ParamCoord(0.75839)}, {2})[0],
              DoubleNear(bspline_1d_after->Evaluate({ParamCoord(0.75839)}, {2})[0], 0.00001));

  ASSERT_THAT(nurbs_1d_->Evaluate({ParamCoord(0.22004)}, {0})[0],
              DoubleEq(nurbs_1d_after->Evaluate({ParamCoord(0.22004)}, {0})[0]));
  ASSERT_THAT(nurbs_1d_->Evaluate({ParamCoord(0.22004)}, {1})[0],
              DoubleEq(nurbs_1d_after->Evaluate({ParamCoord(0.22004)}, {1})[0]));

  ASSERT_THAT(b_spline_2d_->Evaluate({ParamCoord(0.85111)}, {0})[0],
              DoubleEq(bspline_2d_after->Evaluate({ParamCoord(0.85111)}, {0})[0]));
  ASSERT_THAT(b_spline_2d_->Evaluate({ParamCoord(0.85111)}, {1})[0],
              DoubleEq(bspline_2d_after->Evaluate({ParamCoord(0.85111)}, {1})[0]));
  ASSERT_THAT(b_spline_2d_->Evaluate({ParamCoord(0.85111)}, {2})[0],
              DoubleEq(bspline_2d_after->Evaluate({ParamCoord(0.85111)}, {2})[0]));

  ASSERT_THAT(nurbs_2d_->Evaluate({ParamCoord(0.33359)}, {0})[0],
              DoubleEq(nurbs_2d_after->Evaluate({ParamCoord(0.33359)}, {0})[0]));
  ASSERT_THAT(nurbs_2d_->Evaluate({ParamCoord(0.33359)}, {1})[0],
              DoubleEq(nurbs_2d_after->Evaluate({ParamCoord(0.33359)}, {1})[0]));

  ASSERT_THAT(b_spline_3d_->Evaluate({ParamCoord(0.89463)}, {0})[0],
              DoubleEq(bspline_3d_after->Evaluate({ParamCoord(0.89463)}, {0})[0]));
  ASSERT_THAT(b_spline_3d_->Evaluate({ParamCoord(0.89463)}, {1})[0],
              DoubleEq(bspline_3d_after->Evaluate({ParamCoord(0.89463)}, {1})[0]));
  ASSERT_THAT(b_spline_3d_->Evaluate({ParamCoord(0.89463)}, {2})[0],
              DoubleEq(bspline_3d_after->Evaluate({ParamCoord(0.89463)}, {2})[0]));

  ASSERT_THAT(nurbs_3d_->Evaluate({ParamCoord(0.00021)}, {0})[0],
              DoubleEq(nurbs_3d_after->Evaluate({ParamCoord(0.00021)}, {0})[0]));
  ASSERT_THAT(nurbs_3d_->Evaluate({ParamCoord(0.00021)}, {1})[0],
              DoubleEq(nurbs_3d_after->Evaluate({ParamCoord(0.00021)}, {1})[0]));
  ASSERT_THAT(nurbs_3d_->Evaluate({ParamCoord(0.00021)}, {2})[0],
              DoubleEq(nurbs_3d_after->Evaluate({ParamCoord(0.00021)}, {2})[0]));
  remove("splines.itd");
}

TEST_F(AnIRITWriter, ThrowsForSplineOfDimensionFour) {  // NOLINT
  std::shared_ptr<spl::NURBS<4>> nurbs_4d_;
  std::any nurbs_4d_any = std::make_any<std::shared_ptr<spl::NURBS<4>>>(nurbs_4d_);
  ASSERT_THROW(irit_writer_->WriteFile({nurbs_4d_any}, "4d_spline.xml"), std::runtime_error);
  remove("4d_spline.xml");
}
