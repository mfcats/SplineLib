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
    b_spline_ = std::make_shared<spl::BSpline<3>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> b_spline_;
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
    nurbs_ = std::make_shared<spl::NURBS<3>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<3>> nurbs_;
  virtual ~A3DNURBSForIRIT() = default;
};

class A3DIRITReader : public Test, public A3DBSplineForIRIT, public A3DNURBSForIRIT {
 public:
  A3DIRITReader() : irit_reader(std::make_unique<io::IRITReader<3>>()) {}

 protected:
  std::unique_ptr<io::IRITReader<3>> irit_reader;
};

TEST_F(A3DIRITReader, Finds2SplinesOfDimension3) {  // NOLINT
  ASSERT_THAT(irit_reader->ReadIRITFile(path_to_iris_file).size(), 2);
}

TEST_F(A3DIRITReader, ReturnsCorrectDegree) {  // NOLINT
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<3>>>(
      irit_reader->ReadIRITFile(path_to_iris_file)[0])->GetDegree(0).get(), b_spline_->GetDegree(0).get());
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<3>>>(
      irit_reader->ReadIRITFile(path_to_iris_file)[0])->GetDegree(1).get(), b_spline_->GetDegree(1).get());
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<3>>>(
      irit_reader->ReadIRITFile(path_to_iris_file)[0])->GetDegree(2).get(), b_spline_->GetDegree(2).get());

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<3>>>(
      irit_reader->ReadIRITFile(path_to_iris_file)[1])->GetDegree(0).get(), nurbs_->GetDegree(0).get());
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<3>>>(
      irit_reader->ReadIRITFile(path_to_iris_file)[1])->GetDegree(1).get(), nurbs_->GetDegree(1).get());
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<3>>>(
      irit_reader->ReadIRITFile(path_to_iris_file)[1])->GetDegree(2).get(), nurbs_->GetDegree(2).get());
}

TEST_F(A3DIRITReader, ReturnsSameValuesAsGivenSplines) {  // NOLINT
  std::vector<std::any> spline_vector = irit_reader->ReadIRITFile(path_to_iris_file);

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<3>>>(spline_vector[0])->Evaluate({ParamCoord{0.5}}, {0})[0],
              DoubleEq(b_spline_->Evaluate({ParamCoord{0.5}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<3>>>(spline_vector[0])->Evaluate({ParamCoord{0.5}}, {1})[0],
              DoubleEq(b_spline_->Evaluate({ParamCoord{0.5}}, {1})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<3>>>(spline_vector[0])->Evaluate({ParamCoord{0.5}}, {2})[0],
              DoubleEq(b_spline_->Evaluate({ParamCoord{0.5}}, {2})[0]));

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<3>>>(spline_vector[1])->Evaluate({ParamCoord{0.123}}, {0})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.123}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<3>>>(spline_vector[1])->Evaluate({ParamCoord{0.123}}, {1})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.123}}, {1})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<3>>>(spline_vector[1])->Evaluate({ParamCoord{0.123}}, {2})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.123}}, {2})[0]));
}

class A3DIRITWriter : public Test, public A3DBSplineForIRIT, public A3DNURBSForIRIT {
 public:
  A3DIRITWriter() : irit_writer_(std::make_unique<io::IRITWriter<3>>()) {
    std::any b_spline_any = std::make_any<std::shared_ptr<spl::BSpline<3>>>(b_spline_);
    std::any nurbs_any = std::make_any<std::shared_ptr<spl::NURBS<3>>>(nurbs_);
    splines_ = {b_spline_any, nurbs_any};
  }

 protected:
  std::unique_ptr<io::IRITWriter<3>> irit_writer_;
  std::vector<std::any> splines_;
};

TEST_F(A3DIRITWriter, CreatesCorrectFile) {  // NOLINT
  irit_writer_->WriteIRITFile(splines_, "3d_splines.itd");
  std::ifstream newFile;
  newFile.open("3d_splines.itd");
  std::string line, file;
  while (getline(newFile, line)) {
    file += line;
  }
  ASSERT_THAT(file.find("TRIVAR BSPLINE 2 2 2 2 2 2 E3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("TRIVAR BSPLINE 2 2 2 2 2 2 P3"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[KV "), Ne(std::string::npos));
  ASSERT_THAT(file.find(" 1.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[0.000000 1.000000 0.000000]"), Ne(std::string::npos));
  ASSERT_THAT(file.find("[0.750000 1.000000 1.000000 0.000000]"), Ne(std::string::npos));
  remove("3d_splines.itd");
}

TEST_F(A3DIRITWriter, ReturnsSameValuesBeforeAndAfterWritingAndReadingIRITFile) {  // NOLINT
  irit_writer_->WriteIRITFile(splines_, "3d_splines.itd");
  std::unique_ptr<io::IRITReader<3>> irit_reader(std::make_unique<io::IRITReader<3>>());
  auto bspline_after = std::any_cast<std::shared_ptr<spl::BSpline<3>>>(irit_reader->ReadIRITFile("3d_splines.itd")[0]);
  auto nurbs_after = std::any_cast<std::shared_ptr<spl::NURBS<3>>>(irit_reader->ReadIRITFile("3d_splines.itd")[1]);

  ASSERT_THAT(b_spline_->Evaluate({ParamCoord(0.7583), ParamCoord(0.01453), ParamCoord(0.5789)}, {0})[0],
              DoubleEq(bspline_after->Evaluate({ParamCoord(0.7583), ParamCoord(0.01453), ParamCoord(0.5789)}, {0})[0]));
  ASSERT_THAT(b_spline_->Evaluate({ParamCoord(0.7583), ParamCoord(0.01453), ParamCoord(0.5789)}, {1})[0],
              DoubleEq(bspline_after->Evaluate({ParamCoord(0.7583), ParamCoord(0.01453), ParamCoord(0.5789)}, {1})[0]));
  ASSERT_THAT(b_spline_->Evaluate({ParamCoord(0.7583), ParamCoord(0.01453), ParamCoord(0.5789)}, {2})[0],
              DoubleEq(bspline_after->Evaluate({ParamCoord(0.7583), ParamCoord(0.01453), ParamCoord(0.5789)}, {2})[0]));

  ASSERT_THAT(nurbs_->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453), ParamCoord(0.5789)}, {0})[0],
              DoubleEq(nurbs_after->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453), ParamCoord(0.5789)}, {0})[0]));
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453), ParamCoord(0.5789)}, {1})[0],
              DoubleEq(nurbs_after->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453), ParamCoord(0.5789)}, {1})[0]));
  ASSERT_THAT(nurbs_->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453), ParamCoord(0.5789)}, {2})[0],
              DoubleEq(nurbs_after->Evaluate({ParamCoord(0.75839), ParamCoord(0.01453), ParamCoord(0.5789)}, {2})[0]));
  remove("3d_splines.itd");
}

TEST_F(A3DIRITWriter, ReturnsSameValuesBeforeAndAfterConvertingIRITToXMLFile) {  // NOLINT
  io::XMLWriter<3> xml_writer;
  xml_writer.ConvertIRITFileToXMLFile(path_to_iris_file, "converted_xml_file.xml");
  io::XMLReader<3> xml_reader;
  std::vector<std::any> spline_vector = xml_reader.ReadXMLFile("converted_xml_file.xml");
  ASSERT_THAT(spline_vector.size(), 2);

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<3>>>(spline_vector[0])->Evaluate({ParamCoord{0.5}}, {0})[0],
              DoubleEq(b_spline_->Evaluate({ParamCoord{0.5}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<3>>>(spline_vector[0])->Evaluate({ParamCoord{0.5}}, {1})[0],
              DoubleEq(b_spline_->Evaluate({ParamCoord{0.5}}, {1})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::BSpline<3>>>(spline_vector[0])->Evaluate({ParamCoord{0.5}}, {2})[0],
              DoubleEq(b_spline_->Evaluate({ParamCoord{0.5}}, {2})[0]));

  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<3>>>(spline_vector[1])->Evaluate({ParamCoord{0.123}}, {0})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.123}}, {0})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<3>>>(spline_vector[1])->Evaluate({ParamCoord{0.123}}, {1})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.123}}, {1})[0]));
  ASSERT_THAT(std::any_cast<std::shared_ptr<spl::NURBS<3>>>(spline_vector[1])->Evaluate({ParamCoord{0.123}}, {2})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.123}}, {2})[0]));
  remove("converted_xml_file.xml");
}
