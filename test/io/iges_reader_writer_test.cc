/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <config_iges.h>
#include "gmock/gmock.h"
#include "iges_reader.h"
#include "iges_writer.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

using namespace splinelib::src;

class AnIGESReaderAndWriter : public Test {
 public:
  AnIGESReaderAndWriter() {
    baf::KnotVectors<2> knot_vector = {
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.25},
             ParametricCoordinate{0.25},
             ParametricCoordinate{0.5}, ParametricCoordinate{0.5}, ParametricCoordinate{0.75},
             ParametricCoordinate{0.75}, ParametricCoordinate{1},
             ParametricCoordinate{1}, ParametricCoordinate{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.5},
             ParametricCoordinate{0.5},
             ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{1}}))};

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
    b_spline_ = std::make_unique<spl::BSpline<2>>(knot_vector, degree, control_points);
    iges_reader_ = std::make_unique<io::IGESReader>();
    iges_writer_ = std::make_unique<io::IGESWriter>();
  }

 protected:
  std::unique_ptr<spl::NURBS<2>> nurbs_;
  std::unique_ptr<spl::BSpline<2>> b_spline_;
  std::unique_ptr<io::IGESReader> iges_reader_;
  std::unique_ptr<io::IGESWriter> iges_writer_;
};

TEST_F(AnIGESReaderAndWriter, Read1DBSplineFromIGESFile) { // NOLINT
  auto b_spline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_reader_->ReadFile(iges_read)[1]);
  ASSERT_THAT(b_spline_1d->Evaluate({ParametricCoordinate{0.0}}, {0})[0], DoubleNear(-2.23308, 0.0005));
  ASSERT_THAT(b_spline_1d->Evaluate({ParametricCoordinate{0.0}}, {1})[0], DoubleNear(-0.01433, 0.0005));
  ASSERT_THAT(b_spline_1d->Evaluate({ParametricCoordinate{0.0}}, {2})[0], DoubleNear(-0.51255, 0.0005));
  ASSERT_THAT(b_spline_1d->Evaluate({ParametricCoordinate{1.0}}, {0})[0], DoubleNear(-1.3353, 0.0005));
  ASSERT_THAT(b_spline_1d->Evaluate({ParametricCoordinate{1.0}}, {1})[0], DoubleNear(0.450443, 0.0005));
  ASSERT_THAT(b_spline_1d->Evaluate({ParametricCoordinate{1.0}}, {2})[0], DoubleNear(-0.023586, 0.0005));
}

TEST_F(AnIGESReaderAndWriter, Read1DNURBSWithWeigthsOneFromIGESFile) { // NOLINT
  auto nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(iges_reader_->ReadFile(iges_read_2)[1]);
  ASSERT_THAT(nurbs_1d->Evaluate({ParametricCoordinate{0.0}}, {0})[0], DoubleNear(-2.23308, 0.0005));
  ASSERT_THAT(nurbs_1d->Evaluate({ParametricCoordinate{0.0}}, {1})[0], DoubleNear(-0.01433, 0.0005));
  ASSERT_THAT(nurbs_1d->Evaluate({ParametricCoordinate{0.0}}, {2})[0], DoubleNear(-0.51255, 0.0005));
  ASSERT_THAT(nurbs_1d->Evaluate({ParametricCoordinate{1.0}}, {0})[0], DoubleNear(-1.3353, 0.0005));
  ASSERT_THAT(nurbs_1d->Evaluate({ParametricCoordinate{1.0}}, {1})[0], DoubleNear(0.450443, 0.0005));
  ASSERT_THAT(nurbs_1d->Evaluate({ParametricCoordinate{1.0}}, {2})[0], DoubleNear(-0.023586, 0.0005));
}

TEST_F(AnIGESReaderAndWriter, Read2DNURBSFromIGESFile) { // NOLINT
  auto nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_reader_->ReadFile(iges_read)[0]);
  ASSERT_THAT(nurbs_2d->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {0})[0],
              DoubleEq(nurbs_->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {0})[0]));
  ASSERT_THAT(nurbs_2d->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {1})[0],
              DoubleEq(nurbs_->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {1})[0]));
  ASSERT_THAT(nurbs_2d->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {2})[0],
              DoubleEq(nurbs_->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {2})[0]));
  ASSERT_THAT(nurbs_2d->Evaluate({ParametricCoordinate{1.0}, ParametricCoordinate{1.0}}, {0})[0],
              DoubleEq(nurbs_->Evaluate({ParametricCoordinate{1.0}, ParametricCoordinate{1.0}}, {0})[0]));
  ASSERT_THAT(nurbs_2d->Evaluate({ParametricCoordinate{1.0}, ParametricCoordinate{1.0}}, {1})[0],
              DoubleEq(nurbs_->Evaluate({ParametricCoordinate{1.0}, ParametricCoordinate{1.0}}, {1})[0]));
  ASSERT_THAT(nurbs_2d->Evaluate({ParametricCoordinate{1.0}, ParametricCoordinate{1.0}}, {2})[0],
              DoubleEq(nurbs_->Evaluate({ParametricCoordinate{1.0}, ParametricCoordinate{1.0}}, {2})[0]));
}

TEST_F(AnIGESReaderAndWriter, Read2DBSplineFromIGESFile) { // NOLINT
  auto b_spline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_reader_->ReadFile(iges_read_2)[0]);
  ASSERT_THAT(b_spline_2d->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {0})[0],
              DoubleEq(b_spline_->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {0})[0]));
  ASSERT_THAT(b_spline_2d->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {1})[0],
              DoubleEq(b_spline_->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {1})[0]));
  ASSERT_THAT(b_spline_2d->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {2})[0],
              DoubleEq(b_spline_->Evaluate({ParametricCoordinate{0.0}, ParametricCoordinate{0.0}}, {2})[0]));
  ASSERT_THAT(b_spline_2d->Evaluate({ParametricCoordinate{1.0}, ParametricCoordinate{1.0}}, {0})[0],
              DoubleEq(b_spline_->Evaluate({ParametricCoordinate{1.0}, ParametricCoordinate{1.0}}, {0})[0]));
  ASSERT_THAT(b_spline_2d->Evaluate({ParametricCoordinate{1.0}, ParametricCoordinate{1.0}}, {1})[0],
              DoubleEq(b_spline_->Evaluate({ParametricCoordinate{1.0}, ParametricCoordinate{1.0}}, {1})[0]));
  ASSERT_THAT(b_spline_2d->Evaluate({ParametricCoordinate{1.0}, ParametricCoordinate{1.0}}, {2})[0],
              DoubleEq(b_spline_->Evaluate({ParametricCoordinate{1.0}, ParametricCoordinate{1.0}}, {2})[0]));
}

TEST_F(AnIGESReaderAndWriter, Write1DBSplineToIGESFile) { // NOLINT
  auto splines = iges_reader_->ReadFile(iges_read);
  iges_writer_->WriteFile(splines, "write.iges");
  auto b_spline_1d = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(iges_reader_->ReadFile("write.iges")[1]);
  ASSERT_THAT(b_spline_1d->Evaluate({ParametricCoordinate{0.0}}, {0})[0], DoubleNear(-2.23308, 0.0005));
  ASSERT_THAT(b_spline_1d->Evaluate({ParametricCoordinate{0.0}}, {1})[0], DoubleNear(-0.01433, 0.0005));
  ASSERT_THAT(b_spline_1d->Evaluate({ParametricCoordinate{0.0}}, {2})[0], DoubleNear(-0.51255, 0.0005));
  ASSERT_THAT(b_spline_1d->Evaluate({ParametricCoordinate{1.0}}, {0})[0], DoubleNear(-1.3353, 0.0005));
  ASSERT_THAT(b_spline_1d->Evaluate({ParametricCoordinate{1.0}}, {1})[0], DoubleNear(0.450443, 0.0005));
  ASSERT_THAT(b_spline_1d->Evaluate({ParametricCoordinate{1.0}}, {2})[0], DoubleNear(-0.023586, 0.0005));
  remove("write.iges");
}

TEST_F(AnIGESReaderAndWriter, Write1DNURBSToIGESFile) { // NOLINT
  auto splines = iges_reader_->ReadFile(iges_read_2);
  iges_writer_->WriteFile(splines, "write.iges");
  auto nurbs_1d = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(iges_reader_->ReadFile("write.iges")[1]);
  ASSERT_THAT(nurbs_1d->Evaluate({ParametricCoordinate{0.0}}, {0})[0], DoubleNear(-2.23308, 0.0005));
  ASSERT_THAT(nurbs_1d->Evaluate({ParametricCoordinate{0.0}}, {1})[0], DoubleNear(-0.01433, 0.0005));
  ASSERT_THAT(nurbs_1d->Evaluate({ParametricCoordinate{0.0}}, {2})[0], DoubleNear(-0.51255, 0.0005));
  ASSERT_THAT(nurbs_1d->Evaluate({ParametricCoordinate{1.0}}, {0})[0], DoubleNear(-1.3353, 0.0005));
  ASSERT_THAT(nurbs_1d->Evaluate({ParametricCoordinate{1.0}}, {1})[0], DoubleNear(0.450443, 0.0005));
  ASSERT_THAT(nurbs_1d->Evaluate({ParametricCoordinate{1.0}}, {2})[0], DoubleNear(-0.023586, 0.0005));
  remove("write.iges");
}

TEST_F(AnIGESReaderAndWriter, Write2DNURBSToIGESFile) { // NOLINT
  auto splines = iges_reader_->ReadFile(iges_read);
  iges_writer_->WriteFile(splines, "write.iges");
  auto nurbs_2d = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(iges_reader_->ReadFile("write.iges")[0]);
  ASSERT_THAT(nurbs_2d->Evaluate({ParametricCoordinate{0.1}, ParametricCoordinate{0.1}}, {0})[0],
              DoubleNear(nurbs_->Evaluate({ParametricCoordinate{0.1}, ParametricCoordinate{0.1}}, {0})[0], 0.0005));
  ASSERT_THAT(nurbs_2d->Evaluate({ParametricCoordinate{0.1}, ParametricCoordinate{0.1}}, {1})[0],
              DoubleNear(nurbs_->Evaluate({ParametricCoordinate{0.1}, ParametricCoordinate{0.1}}, {1})[0], 0.0005));
  ASSERT_THAT(nurbs_2d->Evaluate({ParametricCoordinate{0.1}, ParametricCoordinate{0.1}}, {2})[0],
              DoubleNear(nurbs_->Evaluate({ParametricCoordinate{0.1}, ParametricCoordinate{0.1}}, {2})[0], 0.0005));
  ASSERT_THAT(nurbs_2d->Evaluate({ParametricCoordinate{0.9}, ParametricCoordinate{0.9}}, {0})[0],
              DoubleNear(nurbs_->Evaluate({ParametricCoordinate{0.9}, ParametricCoordinate{0.9}}, {0})[0], 0.0005));
  ASSERT_THAT(nurbs_2d->Evaluate({ParametricCoordinate{0.9}, ParametricCoordinate{0.9}}, {1})[0],
              DoubleNear(nurbs_->Evaluate({ParametricCoordinate{0.9}, ParametricCoordinate{0.9}}, {1})[0], 0.0005));
  ASSERT_THAT(nurbs_2d->Evaluate({ParametricCoordinate{0.9}, ParametricCoordinate{0.9}}, {2})[0],
              DoubleNear(nurbs_->Evaluate({ParametricCoordinate{0.9}, ParametricCoordinate{0.9}}, {2})[0], 0.0005));
  remove("write.iges");
}

TEST_F(AnIGESReaderAndWriter, Write2DBSplineToIGESFile) { // NOLINT
  auto splines = iges_reader_->ReadFile(iges_read_2);
  iges_writer_->WriteFile(splines, "write.iges");
  auto b_spline_2d = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(iges_reader_->ReadFile("write.iges")[0]);
  ASSERT_THAT(b_spline_2d->Evaluate({ParametricCoordinate{0.1}, ParametricCoordinate{0.1}}, {0})[0],
              DoubleEq(b_spline_->Evaluate({ParametricCoordinate{0.1}, ParametricCoordinate{0.1}}, {0})[0]));
  ASSERT_THAT(b_spline_2d->Evaluate({ParametricCoordinate{0.1}, ParametricCoordinate{0.1}}, {1})[0],
              DoubleEq(b_spline_->Evaluate({ParametricCoordinate{0.1}, ParametricCoordinate{0.1}}, {1})[0]));
  ASSERT_THAT(b_spline_2d->Evaluate({ParametricCoordinate{0.1}, ParametricCoordinate{0.1}}, {2})[0],
              DoubleEq(b_spline_->Evaluate({ParametricCoordinate{0.1}, ParametricCoordinate{0.1}}, {2})[0]));
  ASSERT_THAT(b_spline_2d->Evaluate({ParametricCoordinate{0.9}, ParametricCoordinate{0.9}}, {0})[0],
              DoubleEq(b_spline_->Evaluate({ParametricCoordinate{0.9}, ParametricCoordinate{0.9}}, {0})[0]));
  ASSERT_THAT(b_spline_2d->Evaluate({ParametricCoordinate{0.9}, ParametricCoordinate{0.9}}, {1})[0],
              DoubleNear(b_spline_->Evaluate({ParametricCoordinate{0.9}, ParametricCoordinate{0.9}}, {1})[0], 0.0001));
  ASSERT_THAT(b_spline_2d->Evaluate({ParametricCoordinate{0.9}, ParametricCoordinate{0.9}}, {2})[0],
              DoubleEq(b_spline_->Evaluate({ParametricCoordinate{0.9}, ParametricCoordinate{0.9}}, {2})[0]));
  remove("write.iges");
}

TEST_F(AnIGESReaderAndWriter, ThrowIfFileCantBeOpened) { // NOLINT
  ASSERT_THROW(std::vector<std::any> splines = iges_reader_->ReadFile("a"), std::runtime_error);
}

TEST_F(AnIGESReaderAndWriter, ThrowsForSplineOfDimensionThree) {  // NOLINT
  std::shared_ptr<spl::BSpline<3>> bspline_3d_;
  std::any bspline_3d_any = std::make_any<std::shared_ptr<spl::BSpline<3>>>(bspline_3d_);
  ASSERT_THROW(iges_writer_->WriteFile({bspline_3d_any}, "3d_spline.xml"), std::runtime_error);
  remove("3d_spline.xml");
}
