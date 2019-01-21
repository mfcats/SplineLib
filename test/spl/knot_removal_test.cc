/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <vtk_writer.h>
#include "gmock/gmock.h"

#include "b_spline.h"
#include "nurbs.h"
#include "random_b_spline_generator.h"
#include "random_nurbs_generator.h"
#include "irit_writer.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class BSplineFig5_26 : public Test {  // NOLINT
 public:
  BSplineFig5_26() {
    std::array<Degree, 1> degree = {Degree{3}};
    ParamCoord zero(0), one(1), two(2);
    KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({zero, zero, zero, zero, one, one, one, two, two, two, two}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint({0.0, 0.0}), baf::ControlPoint({0.0, 1.5}), baf::ControlPoint({1.0, 2.0}),
        baf::ControlPoint({2.0, 2.0}), baf::ControlPoint({3.0, 2.0}), baf::ControlPoint({4.0, 1.5}),
        baf::ControlPoint({4.0, 0.0})
    };
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    spl::BSpline<1> b_spline_after(*bspline_1d_before_);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(b_spline_after);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(BSplineFig5_26, RemovesKnot1_0CorrectlyOneTime) {  // NOLINT
  ASSERT_THAT(bspline_1d_after_->RemoveKnot(ParamCoord(1), 0, 0.1), 1);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 1);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(6).get(), DoubleEq(2));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() - 1);
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint({0.0, 0.0}), baf::ControlPoint({0.0, 1.5}), baf::ControlPoint({1.0, 2.0}),
      baf::ControlPoint({3.0, 2.0}), baf::ControlPoint({4.0, 1.5}), baf::ControlPoint({4.0, 0.0})
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, 0), DoubleEq(new_control_points[i].GetValue(0)));
    ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, 1), DoubleEq(new_control_points[i].GetValue(1)));
  }
  for (int i = 0; i <= 50; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(2 * i / 50.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleEq(bspline_1d_before_->Evaluate(param_coord, {0})[0]));
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {1})[0],
                DoubleEq(bspline_1d_before_->Evaluate(param_coord, {1})[0]));
  }
}

TEST_F(BSplineFig5_26, RemovesKnot1_0CorrectlyTwoTimes) {  // NOLINT
  ASSERT_THAT(bspline_1d_after_->RemoveKnot(ParamCoord(1), 0, 0.0, 2), 2);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 2);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(5).get(), DoubleEq(2));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() - 2);
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint({0.0, 0.0}), baf::ControlPoint({0.0, 1.5}), baf::ControlPoint({2.0, 2.5}),
      baf::ControlPoint({4.0, 1.5}), baf::ControlPoint({4.0, 0.0})
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, 0), DoubleEq(new_control_points[i].GetValue(0)));
    ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, 1), DoubleEq(new_control_points[i].GetValue(1)));
  }
  for (int i = 0; i <= 50; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(2 * i / 50.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleEq(bspline_1d_before_->Evaluate(param_coord, {0})[0]));
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {1})[0],
                DoubleEq(bspline_1d_before_->Evaluate(param_coord, {1})[0]));
  }
}

TEST_F(BSplineFig5_26, RemovesKnot1_0CorrectlyThreeTimesAtOnce) {  // NOLINT
  ASSERT_THAT(bspline_1d_after_->RemoveKnot(ParamCoord(1), 0, 0.5, 3), 3);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 3);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(2));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() - 3);
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint({0.0, 0.0}), baf::ControlPoint({0.0, 3.0}), baf::ControlPoint({4.0, 3.0}),
      baf::ControlPoint({4.0, 0.0})
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, 0), DoubleEq(new_control_points[i].GetValue(0)));
    ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, 1), DoubleEq(new_control_points[i].GetValue(1)));
  }
  for (int i = 0; i <= 50; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(2 * i / 50.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleEq(bspline_1d_before_->Evaluate(param_coord, {0})[0]));
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {1})[0],
                DoubleNear(bspline_1d_before_->Evaluate(param_coord, {1})[0], 0.5));
  }
}

TEST_F(BSplineFig5_26, RemovesOnlyTwoKnots1_0CorrectlyWithTolerance0_1) {  // NOLINT
  ASSERT_THAT(bspline_1d_after_->RemoveKnot(ParamCoord(1), 0, 0.1, 3), 2);
}

class BSplineFig5_27 : public Test {  // NOLINT
 public:
  BSplineFig5_27() {
    std::array<Degree, 1> degree = {Degree{3}};
    ParamCoord zero(0), one(1);
    KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({zero, zero, zero, zero, ParamCoord{0.3}, ParamCoord{0.5}, ParamCoord{0.5}, ParamCoord{0.5},
                         ParamCoord{0.7}, ParamCoord{0.7}, one, one, one, one}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint({0.1, 0.0}), baf::ControlPoint({0.0, 1.0}), baf::ControlPoint({1.0, 2.0}),
        baf::ControlPoint({2.5, 2.25}), baf::ControlPoint({3.0, 2.15}), baf::ControlPoint({3.5, 2.0}),
        baf::ControlPoint({4.0, 1.8}), baf::ControlPoint({4.5, 1.0}), baf::ControlPoint({4.3, 0.0}),
        baf::ControlPoint({6.0, 0.0})
    };
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    spl::BSpline<1> b_spline_after(*bspline_1d_before_);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(b_spline_after);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(BSplineFig5_27, RemovesKnot0_3Correctly) {  // NOLINT
  ASSERT_THAT(bspline_1d_after_->RemoveKnot(ParamCoord(0.3), 0, 0.15), 1);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 1);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(0.5));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() - 1);
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint({0.1, 0.0}), baf::ControlPoint({-1.0 / 15.0, 5.0 / 3.0}), baf::ControlPoint({1.75, 2.4}),
      baf::ControlPoint({3.0, 2.15}), baf::ControlPoint({3.5, 2.0}), baf::ControlPoint({4.0, 1.8}),
      baf::ControlPoint({4.5, 1.0}), baf::ControlPoint({4.3, 0.0}), baf::ControlPoint({6.0, 0.0})
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, 0), DoubleEq(new_control_points[i].GetValue(0)));
    ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, 1), DoubleEq(new_control_points[i].GetValue(1)));
  }
  for (int i = 0; i <= 50; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / 50.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleNear(bspline_1d_before_->Evaluate(param_coord, {0})[0], 0.02));
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {1})[0],
                DoubleNear(bspline_1d_before_->Evaluate(param_coord, {1})[0], 0.1));
  }
}

class BSpline2DFig5_28b : public Test {  // NOLINT
 public:
  BSpline2DFig5_28b() {
    std::array<Degree, 2> degree = {Degree{2}, Degree{3}};
    ParamCoord zero(0), one(1);
    KnotVectors<2> knot_vector_before = {
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({zero, zero, zero, ParamCoord{0.25}, ParamCoord{0.5}, ParamCoord{0.75}, one, one, one})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({zero, zero, zero, zero, ParamCoord{0.3}, ParamCoord{0.3}, ParamCoord{0.3}, ParamCoord{0.7},
                             one, one, one, one}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint({0, 1, 4}), baf::ControlPoint({1, 1, 4.5}), baf::ControlPoint({2.5, 1, 3.5}),
        baf::ControlPoint({3.5, 1, 2.5}), baf::ControlPoint({5, 1, 2}), baf::ControlPoint({6.5, 1, 2.5}),

        baf::ControlPoint({0, 1.5, 3.9}), baf::ControlPoint({1, 1.5, 4.4}), baf::ControlPoint({2.5, 1.5, 3.4}),
        baf::ControlPoint({3.5, 1.5, 2.4}), baf::ControlPoint({5, 1.5, 1.9}), baf::ControlPoint({6.5, 1.5, 2.4}),

        baf::ControlPoint({0, 2, 3.8}), baf::ControlPoint({1, 2, 4.3}), baf::ControlPoint({2.5, 2, 3.3}),
        baf::ControlPoint({3.5, 2, 2.3}), baf::ControlPoint({5, 2, 1.8}), baf::ControlPoint({6.5, 2, 2.3}),

        baf::ControlPoint({0, 2.5, 3.7}), baf::ControlPoint({1, 2.5, 4.2}), baf::ControlPoint({2.5, 2.5, 3.2}),
        baf::ControlPoint({3.5, 2.5, 2.2}), baf::ControlPoint({5, 2.5, 1.7}), baf::ControlPoint({6.5, 2.5, 2.2}),

        baf::ControlPoint({0, 3, 3.6}), baf::ControlPoint({1, 3, 4.1}), baf::ControlPoint({2.5, 3, 3.1}),
        baf::ControlPoint({3.5, 3, 2.1}), baf::ControlPoint({5, 3, 1.6}), baf::ControlPoint({6.5, 3, 2.1}),

        baf::ControlPoint({0, 4, 3.3}), baf::ControlPoint({1, 4, 3.8}), baf::ControlPoint({2.5, 4, 2.8}),
        baf::ControlPoint({3.5, 4, 1.8}), baf::ControlPoint({5, 4, 1.3}), baf::ControlPoint({6.5, 4, 1.8}),

        baf::ControlPoint({0, 5, 3}), baf::ControlPoint({1, 5, 3.5}), baf::ControlPoint({2.5, 5, 2.5}),
        baf::ControlPoint({3.5, 5, 1.5}), baf::ControlPoint({5, 5, 1}), baf::ControlPoint({6.5, 5, 1.5}),

        baf::ControlPoint({0, 5.5, 3}), baf::ControlPoint({1, 5.5, 3.5}), baf::ControlPoint({2.5, 5.5, 2.5}),
        baf::ControlPoint({3.5, 5.5, 1.5}), baf::ControlPoint({5, 5.5, 1}), baf::ControlPoint({6.5, 5.5, 1.5})
    };
    bspline_2d_before_ = std::make_shared<spl::BSpline<2>>(knot_vector_before, degree, control_points);
    spl::BSpline<2> b_spline_after(*bspline_2d_before_);
    bspline_2d_after_ = std::make_shared<spl::BSpline<2>>(b_spline_after);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> bspline_2d_before_;
  std::shared_ptr<spl::BSpline<2>> bspline_2d_after_;
};

TEST_F(BSpline2DFig5_28b, RemovesKnot0_3CorrectlyOneTime) {  // NOLINT
  ASSERT_THAT(bspline_2d_after_->RemoveKnot(ParamCoord(0.3), 1, 0.075), 1);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(1)->GetNumberOfKnots() - 1);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(1)->GetKnot(6).get(), DoubleEq(0.7));
  ASSERT_THAT(bspline_2d_after_->GetNumberOfControlPoints(), bspline_2d_before_->GetNumberOfControlPoints() - 6);
  util::MultiIndexHandler<2> coord_handler({31, 31});
  for (int i = 0; i < coord_handler.Get1DLength(); ++i) {
    std::array<ParamCoord, 2> param_coord{ParamCoord(coord_handler[0] / 30.0), ParamCoord(coord_handler[1] / 30.0)};
    for (int j = 0; j < bspline_2d_after_->GetDimension(); ++j) {
      ASSERT_THAT(bspline_2d_after_->Evaluate(param_coord, {j})[0],
                  DoubleNear(bspline_2d_before_->Evaluate(param_coord, {j})[0], 0.075));
    }
  }
}

TEST_F(BSpline2DFig5_28b, RemovesKnot0_3CorrectlyTwoTimes) {  // NOLINT
  ASSERT_THAT(bspline_2d_after_->RemoveKnot(ParamCoord(0.3), 1, 0.12, 2), 2);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(1)->GetNumberOfKnots() - 2);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(1)->GetKnot(5).get(), DoubleEq(0.7));
  ASSERT_THAT(bspline_2d_after_->GetNumberOfControlPoints(), bspline_2d_before_->GetNumberOfControlPoints() - 12);
  util::MultiIndexHandler<2> coord_handler({31, 31});
  for (int i = 0; i < coord_handler.Get1DLength(); ++i) {
    std::array<ParamCoord, 2> param_coord{ParamCoord(coord_handler[0] / 30.0), ParamCoord(coord_handler[1] / 30.0)};
    for (int j = 0; j < bspline_2d_after_->GetDimension(); ++j) {
      ASSERT_THAT(bspline_2d_after_->Evaluate(param_coord, {j})[0],
                  DoubleNear(bspline_2d_before_->Evaluate(param_coord, {j})[0], 0.12));
    }
  }
}

class A3DBSplineToRemoveKnotInDirection2 : public Test {  // NOLINT
 public:
  A3DBSplineToRemoveKnotInDirection2() {
    std::array<Degree, 3> degree = {Degree{2}, Degree{1}, Degree{2}};
    ParamCoord zero(0), one(1);
    KnotVectors<3> knot_vector_before = {
        std::make_shared<baf::KnotVector>(baf::KnotVector({zero, zero, zero, one, one, one})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({zero, zero, one, one})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({zero, zero, zero, ParamCoord{0.3}, ParamCoord{0.3}, one, one, one}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint({0, 1, 4}), baf::ControlPoint({1, 1, 4.5}), baf::ControlPoint({2.5, 1, 3.5}),
        baf::ControlPoint({0, 1, 6}), baf::ControlPoint({1, 1, 6.5}), baf::ControlPoint({2.5, 1, 5.5}),

        baf::ControlPoint({0, 1.5, 3.9}), baf::ControlPoint({1, 1.5, 4.4}), baf::ControlPoint({2.5, 1.5, 3.4}),
        baf::ControlPoint({0, 1.5, 5.9}), baf::ControlPoint({1, 1.5, 6.4}), baf::ControlPoint({2.5, 1.5, 5.4}),

        baf::ControlPoint({0, 2, 3.8}), baf::ControlPoint({1, 2, 4.3}), baf::ControlPoint({2.5, 2, 3.3}),
        baf::ControlPoint({0, 2, 5.8}), baf::ControlPoint({1, 2, 6.3}), baf::ControlPoint({2.5, 2, 5.3}),

        baf::ControlPoint({0, 2.5, 3.7}), baf::ControlPoint({1, 2.5, 4.2}), baf::ControlPoint({2.5, 2.5, 3.2}),
        baf::ControlPoint({0, 2.5, 5.7}), baf::ControlPoint({1, 2.5, 6.2}), baf::ControlPoint({2.5, 2.5, 5.2}),

        baf::ControlPoint({0, 3, 3.6}), baf::ControlPoint({1, 3, 4.1}), baf::ControlPoint({2.5, 3, 3.1}),
        baf::ControlPoint({0, 3, 5.6}), baf::ControlPoint({1, 3, 6.1}), baf::ControlPoint({2.5, 3, 5.1})
    };
    bspline_3d_before_ = std::make_shared<spl::BSpline<3>>(knot_vector_before, degree, control_points);
    spl::BSpline<3> b_spline_after(*bspline_3d_before_);
    bspline_3d_after_ = std::make_shared<spl::BSpline<3>>(b_spline_after);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> bspline_3d_before_;
  std::shared_ptr<spl::BSpline<3>> bspline_3d_after_;
};

TEST_F(A3DBSplineToRemoveKnotInDirection2, RemovesKnot0_3CorrectlyOneTime) {  // NOLINT
  ASSERT_THAT(bspline_3d_after_->RemoveKnot(ParamCoord(0.3), 2, 0.21), 1);
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(2)->GetNumberOfKnots(),
              bspline_3d_before_->GetKnotVector(2)->GetNumberOfKnots() - 1);
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(2)->GetKnot(5).get(), DoubleEq(1));
  ASSERT_THAT(bspline_3d_after_->GetNumberOfControlPoints(), bspline_3d_before_->GetNumberOfControlPoints() - 6);
  util::MultiIndexHandler<3> coord_handler({21, 21, 21});
  for (int i = 0; i < coord_handler.Get1DLength(); ++i) {
    std::array<ParamCoord, 3> param_coord{
        ParamCoord(coord_handler[0] / 20.0), ParamCoord(coord_handler[1] / 20.0), ParamCoord(coord_handler[2] / 20.0)};
    for (int j = 0; j < bspline_3d_after_->GetDimension(); ++j) {
      ASSERT_THAT(bspline_3d_after_->Evaluate(param_coord, {j})[0],
                  DoubleNear(bspline_3d_before_->Evaluate(param_coord, {j})[0], 0.21));
    }
  }
}

TEST_F(A3DBSplineToRemoveKnotInDirection2, RemovesKnot0_3CorrectlyTwoTimes) {  // NOLINT
  ASSERT_THAT(bspline_3d_after_->RemoveKnot(ParamCoord(0.3), 2, 0.21, 2), 2);
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(2)->GetNumberOfKnots(),
              bspline_3d_before_->GetKnotVector(2)->GetNumberOfKnots() - 2);
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(2)->GetKnot(4).get(), DoubleEq(1));
  ASSERT_THAT(bspline_3d_after_->GetNumberOfControlPoints(), bspline_3d_before_->GetNumberOfControlPoints() - 12);
  util::MultiIndexHandler<3> coord_handler({21, 21, 21});
  for (int i = 0; i < coord_handler.Get1DLength(); ++i) {
    std::array<ParamCoord, 3> param_coord{
        ParamCoord(coord_handler[0] / 20.0), ParamCoord(coord_handler[1] / 20.0), ParamCoord(coord_handler[2] / 20.0)};
    for (int j = 0; j < bspline_3d_after_->GetDimension(); ++j) {
      ASSERT_THAT(bspline_3d_after_->Evaluate(param_coord, {j})[0],
                  DoubleNear(bspline_3d_before_->Evaluate(param_coord, {j})[0], 0.281));
    }
  }
}

class NURBSFig5_26 : public Test {  // NOLINT
 public:
  NURBSFig5_26() {
    std::array<Degree, 1> degree = {Degree{3}};
    KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1},
                         ParamCoord{1}, ParamCoord{2}, ParamCoord{2}, ParamCoord{2}, ParamCoord{2}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.5})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.5})),
        baf::ControlPoint(std::vector<double>({4.0, 0.0}))
    };
    std::vector<double> weights = {10, 0.5, 4, 2, 1, 0.5, 4};
    nurbs_1d_before_ = std::make_shared<spl::NURBS<1>>(knot_vector_before, degree, control_points, weights);
    spl::NURBS<1> nurbs_after(*nurbs_1d_before_);
    nurbs_1d_after_ = std::make_shared<spl::NURBS<1>>(nurbs_after);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_before_;
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_after_;
};

TEST_F(NURBSFig5_26, RemovesKnot1_0CorrectlyOneTime) {  // NOLINT
  nurbs_1d_after_->RemoveKnot(ParamCoord(1), 0, 2);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 1);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetKnot(6).get(), DoubleEq(2));
  ASSERT_THAT(nurbs_1d_after_->GetNumberOfControlPoints(), nurbs_1d_before_->GetNumberOfControlPoints() - 1);
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint(std::vector<double>({0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, 1.5})),
      baf::ControlPoint(std::vector<double>({1.0, 2.0})),
      baf::ControlPoint(std::vector<double>({3.0, 2.0})),
      baf::ControlPoint(std::vector<double>({4.0, 1.5})),
      baf::ControlPoint(std::vector<double>({4.0, 0.0}))
  };
  std::vector<double> new_weights = {10, 0.5, 4, 1, 0.5, 4};
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(nurbs_1d_after_->GetControlPoint({i}, j), DoubleEq(new_control_points[i].GetValue(j)));
    }
  }
  double s = 50;
  for (int i = 0; i <= s; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(2 * i / s)};
    ASSERT_THAT(nurbs_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleNear(nurbs_1d_before_->Evaluate(param_coord, {0})[0], 0.7));
    ASSERT_THAT(nurbs_1d_after_->Evaluate(param_coord, {1})[0],
                DoubleNear(nurbs_1d_before_->Evaluate(param_coord, {1})[0], 0.7));
  }
}

class NURBSFig5_27 : public Test {  // NOLINT
 public:
  NURBSFig5_27() {
    std::array<Degree, 1> degree = {Degree{3}};
    KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.5},
                         ParamCoord{0.5}, ParamCoord{0.5}, ParamCoord{0.7}, ParamCoord{0.7}, ParamCoord{1},
                         ParamCoord{1}, ParamCoord{1},
                         ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.1, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.5, 2.25})),
        baf::ControlPoint(std::vector<double>({3.0, 2.15})),
        baf::ControlPoint(std::vector<double>({3.5, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.8})),
        baf::ControlPoint(std::vector<double>({4.5, 1.0})),
        baf::ControlPoint(std::vector<double>({4.3, 0.0})),
        baf::ControlPoint(std::vector<double>({6.0, 0.0}))
    };
    std::vector<double> weights = {1, 0.5, 2, 4, 1, 1, 2, 1, 1, 2};
    nurbs_1d_before_ = std::make_shared<spl::NURBS<1>>(knot_vector_before, degree, control_points, weights);
    spl::NURBS<1> nurbs_after(*nurbs_1d_before_);
    nurbs_1d_after_ = std::make_shared<spl::NURBS<1>>(nurbs_after);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_before_;
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_after_;
};

TEST_F(NURBSFig5_27, RemovesKnot0_5Correctly) {  // NOLINT
  nurbs_1d_after_->RemoveKnot(ParamCoord(0.5), 0, 0.9);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 1);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetKnot(7).get(), DoubleEq(0.7));
  ASSERT_THAT(nurbs_1d_after_->GetNumberOfControlPoints(), nurbs_1d_before_->GetNumberOfControlPoints() - 1);

  io::IRITWriter writer;
  std::any spline1 = std::make_any<std::shared_ptr<spl::NURBS<1>>>(nurbs_1d_before_);
  std::any spline2 = std::make_any<std::shared_ptr<spl::NURBS<1>>>(nurbs_1d_after_);
  writer.WriteFile({spline1, spline2}, "nurbs.itd");
  io::VTKWriter vtkWriter;
  vtkWriter.WriteFile({spline1, spline2}, "nurbs.vtk", {{50}, {50}});
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint(std::vector<double>({0.1, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 2.0})),
      baf::ControlPoint(std::vector<double>({2.5, 124 / 55.0})),
      baf::ControlPoint(std::vector<double>({3.7, 1.96})),
      baf::ControlPoint(std::vector<double>({4.0, 1.8})),
      baf::ControlPoint(std::vector<double>({4.5, 1.0})),
      baf::ControlPoint(std::vector<double>({4.3, 0.0})),
      baf::ControlPoint(std::vector<double>({6.0, 0.0}))
  };
  std::vector<double> new_weights = {1, 0.5, 2, 11.0 / 6.0, 5.0 / 6.0, 2, 1, 1, 2};
  auto a = nurbs_1d_after_->GetControlPoints();
  std::cout << std::endl;
  for (const auto &u : a) {
    std::cout << u << "  ";
  }
  std::cout << std::endl << std::endl;
  auto b = nurbs_1d_after_->GetWeights();
  for (const auto &p : b) {
    std::cout << p << "  ";
  }
  std::cout << std::endl;
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
//      ASSERT_THAT(nurbs_1d_after_->GetControlPoint({i}, j), DoubleEq(new_control_points[i].GetValue(j)));
    }
//    ASSERT_THAT(nurbs_1d_after_->GetWeight({i}), DoubleEq(new_weights[i]));
  }
  double s = 50.0;
  for (int i = 0; i <= s; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / s)};
    ASSERT_THAT(nurbs_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleNear(nurbs_1d_before_->Evaluate(param_coord, {0})[0], 0.9));
    ASSERT_THAT(nurbs_1d_after_->Evaluate(param_coord, {1})[0],
                DoubleNear(nurbs_1d_before_->Evaluate(param_coord, {1})[0], 0.9));
  }
}
