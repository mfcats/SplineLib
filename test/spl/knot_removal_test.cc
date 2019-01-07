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

#include "b_spline.h"
#include "nurbs.h"
#include "random_b_spline_generator.h"
#include "random_nurbs_generator.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class BSplineFig5_26 : public Test {  // NOLINT
 public:
  BSplineFig5_26() {
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
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    spl::BSpline<1> b_spline_after(*bspline_1d_before_);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(b_spline_after);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(BSplineFig5_26, RemovesKnot1_0CorrectlyOneTime) {  // NOLINT
  bspline_1d_after_->RemoveKnot(ParamCoord(1), 0, 0.1);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 1);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(6).get(), DoubleEq(2));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() - 1);
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint(std::vector<double>({0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, 1.5})),
      baf::ControlPoint(std::vector<double>({1.0, 2.0})),
      baf::ControlPoint(std::vector<double>({3.0, 2.0})),
      baf::ControlPoint(std::vector<double>({4.0, 1.5})),
      baf::ControlPoint(std::vector<double>({4.0, 0.0}))
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, j), DoubleEq(new_control_points[i].GetValue(j)));
    }
  }
  double s = 50;
  for (int i = 0; i <= s; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(2 * i / s)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleEq(bspline_1d_before_->Evaluate(param_coord, {0})[0]));
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {1})[0],
                DoubleEq(bspline_1d_before_->Evaluate(param_coord, {1})[0]));
  }
}

TEST_F(BSplineFig5_26, RemovesKnot1_0CorrectlyTwoTimes) {  // NOLINT
  bspline_1d_after_->RemoveKnot(ParamCoord(1), 0, 0.0, 2);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 2);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(5).get(), DoubleEq(2));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() - 2);
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint(std::vector<double>({0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, 1.5})),
      baf::ControlPoint(std::vector<double>({2.0, 2.5})),
      baf::ControlPoint(std::vector<double>({4.0, 1.5})),
      baf::ControlPoint(std::vector<double>({4.0, 0.0}))
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, j), DoubleEq(new_control_points[i].GetValue(j)));
    }
  }
  double s = 50;
  for (int i = 0; i <= s; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(2 * i / s)};
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
      baf::ControlPoint(std::vector<double>({0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, 3.0})),
      baf::ControlPoint(std::vector<double>({4.0, 3.0})),
      baf::ControlPoint(std::vector<double>({4.0, 0.0}))
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, j), DoubleEq(new_control_points[i].GetValue(j)));
    }
  }
  double s = 50;
  for (int i = 0; i <= s; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(2 * i / s)};
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
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    spl::BSpline<1> b_spline_after(*bspline_1d_before_);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(b_spline_after);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(BSplineFig5_27, RemovesKnot0_3Correctly) {  // NOLINT
  bspline_1d_after_->RemoveKnot(ParamCoord(0.3), 0, 0.15);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 1);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(0.5));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() - 1);
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint(std::vector<double>({0.1, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0 / 15.0, 5.0 / 3.0})),
      baf::ControlPoint(std::vector<double>({1.75, 2.4})),
      baf::ControlPoint(std::vector<double>({3.0, 2.15})),
      baf::ControlPoint(std::vector<double>({3.5, 2.0})),
      baf::ControlPoint(std::vector<double>({4.0, 1.8})),
      baf::ControlPoint(std::vector<double>({4.5, 1.0})),
      baf::ControlPoint(std::vector<double>({4.3, 0.0})),
      baf::ControlPoint(std::vector<double>({6.0, 0.0}))
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, j), DoubleEq(new_control_points[i].GetValue(j)));
    }
  }
  double s = 50.0;
  for (int i = 0; i <= s; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / s)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleNear(bspline_1d_before_->Evaluate(param_coord, {0})[0], 0.02));
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {1})[0],
                DoubleNear(bspline_1d_before_->Evaluate(param_coord, {1})[0], 0.1));
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
        baf::ControlPoint(std::vector<double>({0.0, 3.0})),
        baf::ControlPoint(std::vector<double>({0.25, 0.5})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({8.0, 3.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0}))
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
  nurbs_1d_after_->RemoveKnot(ParamCoord(1), 0, 0.1);
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
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(nurbs_1d_after_->GetControlPoint({i}, j), DoubleEq(new_control_points[i].GetValue(j)));
    }
  }
  double s = 50;
  for (int i = 0; i <= s; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(2 * i / s)};
    ASSERT_THAT(nurbs_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleEq(nurbs_1d_before_->Evaluate(param_coord, {0})[0]));
    ASSERT_THAT(nurbs_1d_after_->Evaluate(param_coord, {1})[0],
                DoubleEq(nurbs_1d_before_->Evaluate(param_coord, {1})[0]));
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
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({0.5, 1.0})),
        baf::ControlPoint(std::vector<double>({10.0, 9.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.15})),
        baf::ControlPoint(std::vector<double>({3.5, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.9})),
        baf::ControlPoint(std::vector<double>({4.5, 1.0})),
        baf::ControlPoint(std::vector<double>({4.3, 0.0})),
        baf::ControlPoint(std::vector<double>({3.0, 0.0}))
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

TEST_F(NURBSFig5_27, RemovesKnot0_3Correctly) {  // NOLINT
  nurbs_1d_after_->RemoveKnot(ParamCoord(0.3), 0, 0.15);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 1);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(0.5));
  ASSERT_THAT(nurbs_1d_after_->GetNumberOfControlPoints(), nurbs_1d_before_->GetNumberOfControlPoints() - 1);
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint(std::vector<double>({0.1, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0 / 15.0, 5.0 / 3.0})),
      baf::ControlPoint(std::vector<double>({1.75, 2.4})),
      baf::ControlPoint(std::vector<double>({3.0, 2.15})),
      baf::ControlPoint(std::vector<double>({3.5, 2.0})),
      baf::ControlPoint(std::vector<double>({4.0, 1.8})),
      baf::ControlPoint(std::vector<double>({4.5, 1.0})),
      baf::ControlPoint(std::vector<double>({4.3, 0.0})),
      baf::ControlPoint(std::vector<double>({6.0, 0.0}))
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(nurbs_1d_after_->GetControlPoint({i}, j), DoubleEq(new_control_points[i].GetValue(j)));
    }
  }
  double s = 50.0;
  for (int i = 0; i <= s; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / s)};
    ASSERT_THAT(nurbs_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleNear(nurbs_1d_before_->Evaluate(param_coord, {0})[0], 0.02));
    ASSERT_THAT(nurbs_1d_after_->Evaluate(param_coord, {1})[0],
                DoubleNear(nurbs_1d_before_->Evaluate(param_coord, {1})[0], 0.1));
  }
}
