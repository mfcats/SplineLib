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

class A1DBNURBSForVTK : public Test {  // NOLINT
 public:
  A1DBNURBSForVTK() : vtk_writer_(std::make_unique<io::VTKWriter>()) {
    std::array<Degree, 1> degree = {Degree{1}};
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0})),
        baf::ControlPoint(std::vector<double>({1.0})),
        baf::ControlPoint(std::vector<double>({1.2}))
    };
    std::vector<double> weights = {3.0, 0.5, 1.8};
    nurbs_1d_ = std::make_shared<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_;
  std::unique_ptr<io::VTKWriter> vtk_writer_;
};

TEST_F(A1DBNURBSForVTK, CreatesVTKFile) {
  std::any nurbs_1d_any = std::make_any<std::shared_ptr<spl::NURBS<1>>>(nurbs_1d_);
  vtk_writer_->WriteFile({nurbs_1d_any}, "test.vtk", {{10}});
  std::ifstream newFile;
  newFile.open("test.vtk");
  if (!newFile.good()) {
    throw std::runtime_error("VTK file could not be opened.");
  }
  std::string line;
  std::cout << std::endl << std::endl;
  while (getline(newFile, line)) {
    std::cout << line << std::endl;
  }
  std::cout << std::endl;
  // remove("test.vtk");
}