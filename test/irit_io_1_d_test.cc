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
    b_spline_1_ = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> b_spline_1_;
  virtual ~A1DBSplineForIRIT() = default;
};

class ASecond1DBSplineForIRIT {  // NOLINT
 public:
  ASecond1DBSplineForIRIT() {
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
    b_spline_2_ = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> b_spline_2_;
  virtual ~ASecond1DBSplineForIRIT() = default;
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
    nurbs_ = std::make_shared<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_;
  virtual ~A1DNURBSForIRIT() = default;
};

class A1DIRITReader : public Test, public A1DBSplineForIRIT, public ASecond1DBSplineForIRIT, public A1DNURBSForIRIT {
 public:
  A1DIRITReader() : irit_reader(std::make_unique<io::IRITReader<1>>()) {}

 protected:
  std::unique_ptr<io::IRITReader<1>> irit_reader;
};

TEST_F(A1DIRITReader, ThrowsExceptionForNonExistingFile) {  // NOLINT
  ASSERT_THROW(irit_reader->ReadIRITFile("testing.itd"), std::runtime_error);
}

TEST_F(A1DIRITReader, Finds3SplinesOfDimension1) {  // NOLINT
  ASSERT_THAT(irit_reader->ReadIRITFile(path_to_iris_file).size(), 3);
}

TEST_F(A1DIRITReader, ReturnsCorrectDegree) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<1>>>(
      irit_reader->ReadIRITFile(path_to_iris_file)[0])->GetDegree(0).get(), b_spline_1_->GetDegree(0).get());
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<1>>>(
      irit_reader->ReadIRITFile(path_to_iris_file)[1])->GetDegree(0).get(), b_spline_2_->GetDegree(0).get());
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<1>>>(
      irit_reader->ReadIRITFile(path_to_iris_file)[2])->GetDegree(0).get(), nurbs_->GetDegree(0).get());
}

TEST_F(A1DIRITReader, ReturnsSameValuesAsGivenSplines) {  // NOLINT
  std::vector<std::any> spline_vector = irit_reader->ReadIRITFile(path_to_iris_file);
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<1>>>(spline_vector[0])->Evaluate({ParamCoord{10.5}}, {0})[0],
              DoubleEq(b_spline_1_->Evaluate({ParamCoord{10.5}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<1>>>(spline_vector[1])->Evaluate({ParamCoord{0.5}}, {0})[0],
              DoubleEq(b_spline_2_->Evaluate({ParamCoord{0.5}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<1>>>(spline_vector[2])->Evaluate({ParamCoord{0.123}}, {0})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.123}}, {0})[0]));
}

class A1DIRITWriter : public Test, public A1DBSplineForIRIT, public ASecond1DBSplineForIRIT, public A1DNURBSForIRIT {
 public:
  A1DIRITWriter() : irit_writer_(std::make_unique<io::IRITWriter<1>>()) {
    std::any b_spline_1_any = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_1_);
    std::any b_spline_2_any = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_2_);
    std::any nurbs_any = std::make_any<std::shared_ptr<spl::NURBS<1>>>(nurbs_);
    splines_ = {b_spline_1_any, b_spline_2_any, nurbs_any};
  }

 protected:
  std::unique_ptr<io::IRITWriter<1>> irit_writer_;
  std::vector<std::any> splines_;
};

TEST_F(A1DIRITWriter, ThrowsExceptionForWrongInputType) {  // NOLINT
  ASSERT_THROW(io::IRITWriter<1>().WriteIRITFile({std::make_any<int>(8)}, "1d_splines.itd"), std::runtime_error);
  remove("1d_splines.itd");
}

TEST_F(A1DIRITWriter, CreatesCorrectFile) {  // NOLINT
  irit_writer_->WriteIRITFile(splines_, "1d_splines.itd");
  std::ifstream newFile;
  newFile.open("1d_splines.itd");
  std::string line, file;
  while (getline(newFile, line)) {
    file += line;
  }
  ASSERT_THAT(file.find("CURVE BSPLINE 16 4 E2"), Ne(std::string::npos));
  ASSERT_THAT(file.find("CURVE BSPLINE 8 3 E3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("CURVE BSPLINE 6 3 P2"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[KV "), Ne(std::string::npos));
  ASSERT_THAT(file.find(" 11.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find(" 1.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[0.000000 3.572000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[0.800000 0.300000 -0.400000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[0.500000 2.000000 0.500000]"), Ne(std::string::npos));
  remove("1d_splines.itd");
}

TEST_F(A1DIRITWriter, ReturnsSameValuesBeforeAndAfterWritingAndReadingIRITFile) {  // NOLINT
  irit_writer_->WriteIRITFile(splines_, "1d_splines.itd");
  std::unique_ptr<io::IRITReader<1>> irit_reader(std::make_unique<io::IRITReader<1>>());
  auto bspline1_after = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(irit_reader->ReadIRITFile("1d_splines.itd")[0]);
  auto bspline2_after = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(irit_reader->ReadIRITFile("1d_splines.itd")[1]);
  auto nurbs_after = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(irit_reader->ReadIRITFile("1d_splines.itd")[2]);
  ASSERT_THAT(b_spline_1_->Evaluate({ParamCoord(0.75839)}, {0})[0],
              DoubleEq(bspline1_after->Evaluate({ParamCoord(0.75839)}, {0})[0]));
  ASSERT_THAT(b_spline_2_->Evaluate({ParamCoord(0.17456)}, {0})[0],
              DoubleNear(bspline2_after->Evaluate({ParamCoord(0.17456)}, {0})[0], 0.000001));
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord(0.48752)}, {0})[0],
              DoubleEq(nurbs_after->Evaluate({ParamCoord(0.48752)}, {0})[0]));
  remove("1d_splines.itd");
}

TEST_F(A1DIRITWriter, ReturnsSameValuesBeforeAndAfterConvertingIRITToXMLFile) {  // NOLINT
  io::XMLWriter<1> xml_writer;
  xml_writer.ConvertIRITFileToXMLFile(path_to_iris_file, "converted_xml_file.xml");
  io::XMLReader<1> xml_reader;
  std::vector<std::any> spline_vector = xml_reader.ReadXMLFile("converted_xml_file.xml");
  ASSERT_THAT(spline_vector.size(), 3);

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<1>>>(spline_vector[0])->Evaluate({ParamCoord{10.5}}, {0})[0],
              DoubleEq(b_spline_1_->Evaluate({ParamCoord{10.5}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<1>>>(spline_vector[1])->Evaluate({ParamCoord{0.5}}, {0})[0],
              DoubleEq(b_spline_2_->Evaluate({ParamCoord{0.5}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<1>>>(spline_vector[2])->Evaluate({ParamCoord{0.123}}, {0})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.123}}, {0})[0]));
  remove("converted_xml_file.xml");
}
