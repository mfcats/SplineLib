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
using testing::Ne;

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
    b_spline_ = std::make_shared<spl::BSpline<2>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> b_spline_;
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
    nurbs_ = std::make_shared<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<2>> nurbs_;
  virtual ~A2DNURBSForIRIT() = default;
};

class A2DIRITReader : public Test, public A2DBSplineForIRIT, public A2DNURBSForIRIT {
 public:
  A2DIRITReader() : irit_reader(std::make_unique<io::IRITReader<2>>()) {}

 protected:
  std::unique_ptr<io::IRITReader<2>> irit_reader;
};

TEST_F(A2DIRITReader, Finds2SplinesOfDimension2) {  // NOLINT
  ASSERT_THAT(irit_reader->ReadIRITFile(path_to_iris_file).size(), 2);
}

TEST_F(A2DIRITReader, ReturnsCorrectDegree) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      irit_reader->ReadIRITFile(path_to_iris_file)[0])->GetDegree(0).get(), b_spline_->GetDegree(0).get());
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(
      irit_reader->ReadIRITFile(path_to_iris_file)[0])->GetDegree(1).get(), b_spline_->GetDegree(1).get());

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      irit_reader->ReadIRITFile(path_to_iris_file)[1])->GetDegree(0).get(), nurbs_->GetDegree(0).get());
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(
      irit_reader->ReadIRITFile(path_to_iris_file)[1])->GetDegree(1).get(), nurbs_->GetDegree(1).get());
}

TEST_F(A2DIRITReader, ReturnsSameValuesAsGivenSplines) {  // NOLINT
  std::vector<std::any> spline_vector = irit_reader->ReadIRITFile(path_to_iris_file);

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(spline_vector[0])->Evaluate({ParamCoord{0.5}}, {0})[0],
              DoubleEq(b_spline_->Evaluate({ParamCoord{0.5}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(spline_vector[0])->Evaluate({ParamCoord{0.5}}, {1})[0],
              DoubleEq(b_spline_->Evaluate({ParamCoord{0.5}}, {1})[0]));

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(spline_vector[1])->Evaluate({ParamCoord{0.123}}, {0})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.123}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(spline_vector[1])->Evaluate({ParamCoord{0.123}}, {1})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.123}}, {1})[0]));
}

class A2DIRITWriter : public Test, public A2DBSplineForIRIT, public A2DNURBSForIRIT {
 public:
  A2DIRITWriter() : irit_writer_(std::make_unique<io::IRITWriter<2>>()) {
    std::any b_spline_any = std::make_any<std::shared_ptr<spl::BSpline<2>>>(b_spline_);
    std::any nurbs_any = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbs_);
    splines_ = {b_spline_any, nurbs_any};
  }

 protected:
  std::unique_ptr<io::IRITWriter<2>> irit_writer_;
  std::vector<std::any> splines_;
};

TEST_F(A2DIRITWriter, CreatesCorrectFile) {  // NOLINT
  irit_writer_->WriteIRITFile(splines_, "2d_splines.itd");
  std::ifstream newFile;
  newFile.open("2d_splines.itd");
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
  remove("2d_splines.itd");
}

TEST_F(A2DIRITWriter, ReturnsSameValuesBeforeAndAfterWritingAndReadingIRITFile) {  // NOLINT
  irit_writer_->WriteIRITFile(splines_, "2d_splines.itd");
  std::unique_ptr<io::IRITReader<2>> irit_reader(std::make_unique<io::IRITReader<2>>());
  auto bspline_after = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(irit_reader->ReadIRITFile("2d_splines.itd")[0]);
  auto nurbs_after = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(irit_reader->ReadIRITFile("2d_splines.itd")[1]);

  ASSERT_THAT(b_spline_->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453)}, {0})[0],
              DoubleEq(bspline_after->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453)}, {0})[0]));
  ASSERT_THAT(b_spline_->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453)}, {1})[0],
              DoubleEq(bspline_after->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453)}, {1})[0]));

  ASSERT_THAT(nurbs_->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453)}, {0})[0],
              DoubleEq(nurbs_after->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453)}, {0})[0]));
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453)}, {1})[0],
              DoubleEq(nurbs_after->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453)}, {1})[0]));
  remove("2d_splines.itd");
}

TEST_F(A2DIRITWriter, ReturnsSameValuesBeforeAndAfterConvertingIRITToXMLFile) {  // NOLINT
  io::XMLWriter<2> xml_writer;
  xml_writer.ConvertIRITFileToXMLFile(path_to_iris_file, "converted_xml_file.xml");
  io::XMLReader<2> xml_reader;
  std::vector<std::any> spline_vector = xml_reader.ReadXMLFile("converted_xml_file.xml");
  ASSERT_THAT(spline_vector.size(), 2);

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(spline_vector[0])->Evaluate({ParamCoord{0.5}}, {0})[0],
              DoubleEq(b_spline_->Evaluate({ParamCoord{0.5}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<2>>>(spline_vector[0])->Evaluate({ParamCoord{0.5}}, {1})[0],
              DoubleEq(b_spline_->Evaluate({ParamCoord{0.5}}, {1})[0]));

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(spline_vector[1])->Evaluate({ParamCoord{0.123}}, {0})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.123}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<2>>>(spline_vector[1])->Evaluate({ParamCoord{0.123}}, {1})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.123}}, {1})[0]));
  remove("converted_xml_file.xml");
}
