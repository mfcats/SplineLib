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

#include "vtk_writer.h"

using testing::Test;
using testing::Ne;

class A1DBNURBSForVTKWriter {  // NOLINT
 public:
  A1DBNURBSForVTKWriter() {
    std::array<Degree, 1> degree = {Degree{1}};
    KnotVectors<1> knot_vector = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({4.0, -1.0})),
        baf::ControlPoint(std::vector<double>({5.0, 0.0})),
        baf::ControlPoint(std::vector<double>({5.2, 2.0}))
    };
    std::vector<double> weights = {3.0, 0.5, 1.8};
    nurbs_1d_ = std::make_shared<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_;
};

class A2DBNURBSForVTKWriter {  // NOLINT
 public:
  A2DBNURBSForVTKWriter() {
    std::array<Degree, 2> degree = {Degree{1}, Degree{1}};
    KnotVectors<2> knot_vector = {
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({-5.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({-4.0, 0.3, 2.0})),
        baf::ControlPoint(std::vector<double>({-3.8, 0.2, 1.5})),
        baf::ControlPoint(std::vector<double>({-5.5, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({-3.7, 0.8, -0.5})),
        baf::ControlPoint(std::vector<double>({-3.5, 0.9, 0.5})),
        baf::ControlPoint(std::vector<double>({-4.9, 1.2, -1.0})),
        baf::ControlPoint(std::vector<double>({-4.2, 1.5, -0.5})),
        baf::ControlPoint(std::vector<double>({-3.5, 2.0, -1.5}))
    };
    std::vector<double> weights = {3.0, 0.5, 1.8, 2.0, 2.5, 1.8, 1.0, 0.8, 0.2};
    nurbs_2d_ = std::make_shared<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<2>> nurbs_2d_;
};

class A3DBNURBSForVTKWriter {  // NOLINT
 public:
  A3DBNURBSForVTKWriter() {
    std::array<Degree, 3> degree = {Degree{1}, Degree{1}, Degree{1}};
    KnotVectors<3> knot_vector = {
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.3, -0.1})),
        baf::ControlPoint(std::vector<double>({0.2, 2.2, 0.2})),
        baf::ControlPoint(std::vector<double>({1.5, 2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.3, -0.1, 1.2})),
        baf::ControlPoint(std::vector<double>({1.5, 0.1, 1.8})),
        baf::ControlPoint(std::vector<double>({0.1, 2.2, 1.6})),
        baf::ControlPoint(std::vector<double>({1.7, 2.0, 2.2}))
    };
    std::vector<double> weights = {3.0, 0.5, 1.8, 2.0, 1.8, 1.0, 0.8, 0.2};
    nurbs_3d_ = std::make_shared<spl::NURBS<3>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<3>> nurbs_3d_;
};

class AVTKWriter : public Test,
                   public A1DBNURBSForVTKWriter, public A2DBNURBSForVTKWriter, public A3DBNURBSForVTKWriter {  // NOLINT
 public:
  AVTKWriter() : vtk_writer_(std::make_unique<io::VTKWriter>()), scattering_({{100}, {30, 30}, {10, 8, 5}}) {
    std::any nurbs_1d_any = std::make_any<std::shared_ptr<spl::NURBS<1>>>(nurbs_1d_);
    std::any nurbs_2d_any = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbs_2d_);
    std::any nurbs_3d_any = std::make_any<std::shared_ptr<spl::NURBS<3>>>(nurbs_3d_);
    splines_ = {nurbs_1d_any, nurbs_2d_any, nurbs_3d_any};
  }

 protected:
  std::unique_ptr<io::VTKWriter> vtk_writer_;
  std::vector<std::any> splines_;
  std::vector<std::vector<int>> scattering_;
};

TEST_F(AVTKWriter, CreatesVTKFile) {  // NOLINT
  vtk_writer_->WriteFile(splines_, "splines.vtk", scattering_);
  std::ifstream newFile;
  newFile.open("splines.vtk");
  ASSERT_THAT(newFile.good(), true);
  std::string line, file;
  while (getline(newFile, line)) {
    file += line + "\n";
  }
  ASSERT_THAT(file.find("# vtk DataFile Version 3.0\nSpline from Splinelib\nASCII\n"), Ne(std::string::npos));
  ASSERT_THAT(file.find("DATASET UNSTRUCTURED_GRID\nPOINTS 1656 double\n"), Ne(std::string::npos));
  ASSERT_THAT(file.find("CELLS 1400 8400\n"), Ne(std::string::npos));
  ASSERT_THAT(file.find("CELL_TYPES 1400\n"), Ne(std::string::npos));
  remove("splines.vtk");
}

TEST_F(AVTKWriter, ThrowsForMissingEntryInScattering) {  // NOLINT
  scattering_.pop_back();
  ASSERT_THROW(vtk_writer_->WriteFile(splines_, "splines.vtk", scattering_), std::runtime_error);
}

TEST_F(AVTKWriter, ThrowsForTooManyEntriesInScattering) {  // NOLINT
  scattering_.push_back({10});
  ASSERT_THROW(vtk_writer_->WriteFile(splines_, "splines.vtk", scattering_), std::runtime_error);
}

TEST_F(AVTKWriter, ThrowsForMissingEntryInSplineScattering) {  // NOLINT
  scattering_.back().pop_back();
  ASSERT_THROW(vtk_writer_->WriteFile(splines_, "splines.vtk", scattering_), std::runtime_error);
}

TEST_F(AVTKWriter, ThrowsForTooManyEntriesInSplineScattering) {  // NOLINT
  scattering_[1].push_back(20);
  ASSERT_THROW(vtk_writer_->WriteFile(splines_, "splines.vtk", scattering_), std::runtime_error);
}

TEST_F(AVTKWriter, ThrowsForSplineOfDimensionFour) {  // NOLINT
  std::shared_ptr<spl::NURBS<4>> nurbs_4d_;
  std::any nurbs_4d_any = std::make_any<std::shared_ptr<spl::NURBS<4>>>(nurbs_4d_);
  ASSERT_THROW(vtk_writer_->WriteFile({nurbs_4d_any}, "4d_spline.vtk", {{10, 10, 10, 10}}), std::runtime_error);
  remove("4d_spline.vtk");
}
