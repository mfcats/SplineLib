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

class A1DBSplineFromIRITFile : public Test {
 public:
  A1DBSplineFromIRITFile() : irit_reader(std::make_unique<io::IRITReader<1>>()) {
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
        baf::ControlPoint(std::vector<double>({0, 3.572})),
    };

    b_spline_ = std::make_unique<spl::BSpline<1>>(knot_vector, degree, control_points);
  }

 protected:
  std::unique_ptr<spl::BSpline<1>> b_spline_;
  std::unique_ptr<io::IRITReader<1>> irit_reader;
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
