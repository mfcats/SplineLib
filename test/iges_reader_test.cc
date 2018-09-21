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
#include "iges_reader.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class A2DNurbsFromIGESFile : public Test {
 public:
  A2DNurbsFromIGESFile() {
    std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector =
        {std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                            ParamCoord{0.25},
                                                            ParamCoord{0.25}, ParamCoord{0.5}, ParamCoord{0.5},
                                                            ParamCoord{0.75},
                                                            ParamCoord{0.75}, ParamCoord{1}, ParamCoord{1},
                                                            ParamCoord{1}})),
         std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                            ParamCoord{0.5},
                                                            ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1},
                                                            ParamCoord{1}}))};

    std::array<Degree, 2> degree = {Degree{2}, Degree{2}};

    std::vector<double> weights = {1.0, 0.70710678118654757, 1.0, 0.70710678118654757, 1.0, 0.70710678118654757,
                                   1.0, 0.70710678118654757, 1.0, 0.70710678118654757, 0.50000000000000011,
                                   0.70710678118654757, 0.50000000000000011, 0.70710678118654757, 0.50000000000000011,
                                   0.70710678118654757, 0.50000000000000011, 0.70710678118654757, 1.0,
                                   0.70710678118654757, 1.0, 0.70710678118654757, 1.0, 0.70710678118654757, 1.0,
                                   0.70710678118654757, 1.0, 0.70710678118654757, 0.50000000000000011,
                                   0.70710678118654757, 0.50000000000000011, 0.70710678118654757, 0.50000000000000011,
                                   0.70710678118654757, 0.50000000000000011, 0.70710678118654757, 1.0,
                                   0.70710678118654757, 1.0, 0.70710678118654757, 1.0, 0.70710678118654757, 1.0,
                                   0.70710678118654757, 1.0};

    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0, -1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, -0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, -1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, -0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, -1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, -0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, -1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, -0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, -1.0})),
        baf::ControlPoint(std::vector<double>({0.99999999999999989, 0.0, -0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({0.99999999999999978, 0.99999999999999944, -0.99999999999999956})),
        baf::ControlPoint(std::vector<double>({6.1232339957367648E-17, 0.99999999999999989, -0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({-0.99999999999999956, 0.99999999999999978, -0.99999999999999956})),
        baf::ControlPoint(std::vector<double>({-0.99999999999999989, 1.224646799147353E-16, -0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({-0.99999999999999978, -0.99999999999999956, -0.99999999999999956})),
        baf::ControlPoint(std::vector<double>({-1.8369701987210294E-16, -0.99999999999999989, -0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({0.99999999999999978, -0.99999999999999956, -0.99999999999999956})),
        baf::ControlPoint(std::vector<double>({0.99999999999999989, 0.0, -0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.99999999999999989, 0.99999999999999989, 0.0})),
        baf::ControlPoint(std::vector<double>({6.123233995736766E-17, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({-0.99999999999999989, 0.99999999999999967, 0.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 1.2246467991473532E-16, 0.0})),
        baf::ControlPoint(std::vector<double>({-1.0, -0.99999999999999989, 0.0})),
        baf::ControlPoint(std::vector<double>({-1.8369701987210297E-16, -1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, -1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.99999999999999989, 0.0, 0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({0.99999999999999978, 0.99999999999999944, 0.99999999999999956})),
        baf::ControlPoint(std::vector<double>({6.1232339957367648E-17, 0.99999999999999989, 0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({-0.99999999999999956, 0.99999999999999978, 0.99999999999999956})),
        baf::ControlPoint(std::vector<double>({-0.99999999999999989, 1.224646799147353E-16, 0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({-0.99999999999999978, -0.99999999999999956, 0.99999999999999956})),
        baf::ControlPoint(std::vector<double>({-1.8369701987210294E-16, -0.99999999999999989, 0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({0.99999999999999978, -0.99999999999999956, 0.99999999999999956})),
        baf::ControlPoint(std::vector<double>({0.99999999999999989, 0.0, 0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.99999999999999989})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 1.0}))
    };

    nurbs_ = std::make_unique<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
    iges_reader_ = std::make_unique<io::IGESReader>();
  }

 protected:
  std::unique_ptr<spl::NURBS<2>> nurbs_;
  std::unique_ptr<io::IGESReader> iges_reader_;
};

TEST_F(A2DNurbsFromIGESFile, Read1DBSplineFromIGESFile) { // NOLINT
  spl::BSpline<1> b_spline = std::any_cast<spl::BSpline<1>>(iges_reader_->ReadIGESFile(path_to_iges_file)[1]);
  ASSERT_THAT(b_spline.Evaluate({ParamCoord{0.0}}, {0})[0], DoubleNear(-2.23308, 0.0005));
  ASSERT_THAT(b_spline.Evaluate({ParamCoord{0.0}}, {1})[0], DoubleNear(-0.01433, 0.0005));
  ASSERT_THAT(b_spline.Evaluate({ParamCoord{0.0}}, {2})[0], DoubleNear(-0.51255, 0.0005));
  ASSERT_THAT(b_spline.Evaluate({ParamCoord{1.0}}, {0})[0], DoubleNear(-1.3353, 0.0005));
  ASSERT_THAT(b_spline.Evaluate({ParamCoord{1.0}}, {1})[0], DoubleNear(0.450443, 0.0005));
  ASSERT_THAT(b_spline.Evaluate({ParamCoord{1.0}}, {2})[0], DoubleNear(-0.023586, 0.0005));
}

TEST_F(A2DNurbsFromIGESFile, Read2DNURBSFromIGESFile) { // NOLINT
  spl::NURBS<2> nurbs = std::any_cast<spl::NURBS<2>>(iges_reader_->ReadIGESFile(path_to_iges_file)[0]);
  ASSERT_THAT(nurbs.Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {0})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {0})[0]));
  ASSERT_THAT(nurbs.Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {1})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {1})[0]));
  ASSERT_THAT(nurbs.Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {2})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {2})[0]));
  ASSERT_THAT(nurbs.Evaluate({ParamCoord{1.0}, ParamCoord{1.0}}, {0})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{1.0}, ParamCoord{1.0}}, {0})[0]));
  ASSERT_THAT(nurbs.Evaluate({ParamCoord{1.0}, ParamCoord{1.0}}, {1})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{1.0}, ParamCoord{1.0}}, {1})[0]));
  ASSERT_THAT(nurbs.Evaluate({ParamCoord{1.0}, ParamCoord{1.0}}, {2})[0],
              DoubleEq(nurbs_->Evaluate({ParamCoord{1.0}, ParamCoord{1.0}}, {2})[0]));
}

TEST_F(A2DNurbsFromIGESFile, ThrowIfFileCantBeOpened) { // NOLINT
  std::vector<std::any> splines = iges_reader_->ReadIGESFile("a");
}
