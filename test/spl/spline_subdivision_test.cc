/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <any>
#include "gmock/gmock.h"

#include "b_spline.h"
#include "nurbs.h"
#include "random_b_spline_generator.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class BSplineFig5_9 : public Test {  // NOLINT
 public:
  BSplineFig5_9() {
    std::array<Degree, 2> degree = {Degree{3}, Degree{2}};
    std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector = {
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                           ParamCoord{1}, ParamCoord{1}, ParamCoord{1},
                                                           ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5},
                                                           ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({5.0, 0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 3.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 0.0, 3.0})),

        baf::ControlPoint(std::vector<double>({5.0, 2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.5, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 3.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 1.0, 3.0})),

        baf::ControlPoint(std::vector<double>({5.0, 2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.5, 2.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 2.5, 2.0})),

        baf::ControlPoint(std::vector<double>({5.0, 5.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 4.5, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 4.0, 2.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 4.0, 2.0}))
    };
    bspline_2d_ = std::make_shared<spl::BSpline<2>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> bspline_2d_;
};

TEST_F(BSplineFig5_9, IsSubdividedAtKnot0_3InFirstDirection) {  // NOLINT
  auto splines = bspline_2d_->SudivideSpline(ParamCoord{0.7}, 0);
  auto spline1 = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(splines[0]);
  auto spline2 = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(splines[1]);
  ASSERT_THAT(spline1->GetKnotVector(0)->GetNumberOfKnots(), 8);
  ASSERT_THAT(spline2->GetKnotVector(0)->GetNumberOfKnots(), 8);
  ASSERT_THAT(spline1->GetKnotVector(0)->GetKnot(0).get(), DoubleEq(0));
  ASSERT_THAT(spline1->GetKnotVector(0)->GetKnot(3).get(), DoubleEq(0));
  ASSERT_THAT(spline1->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(0.7));
  ASSERT_THAT(spline1->GetKnotVector(0)->GetKnot(7).get(), DoubleEq(0.7));
  ASSERT_THAT(spline2->GetKnotVector(0)->GetKnot(0).get(), DoubleEq(0.7));
  ASSERT_THAT(spline2->GetKnotVector(0)->GetKnot(3).get(), DoubleEq(0.7));
  ASSERT_THAT(spline2->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(1));
  ASSERT_THAT(spline2->GetKnotVector(0)->GetKnot(7).get(), DoubleEq(1));
  for (int i = 0; i <= 100; ++i) {
    for (int j = 0; j <= 100; ++j) {
      std::array<ParamCoord, 2> param_coord{ParamCoord(i / 100.0), ParamCoord(j / 100.0)};
      if (i / 100.0 <= 0.7) {
        ASSERT_THAT(spline1->Evaluate(param_coord, {0})[0],
                    DoubleNear(bspline_2d_->Evaluate(param_coord, {0})[0], 0.000001));
      }
      if (i / 100.0 >= 0.7) {
        ASSERT_THAT(spline2->Evaluate(param_coord, {0})[0],
                    DoubleNear(bspline_2d_->Evaluate(param_coord, {0})[0], 0.000001));
      }
    }
  }
}

TEST_F(BSplineFig5_9, IsSubdividedAtKnot0_4InSecondDirection) {  // NOLINT
  auto splines = bspline_2d_->SudivideSpline(ParamCoord{0.4}, 1);
  auto spline1 = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(splines[0]);
  auto spline2 = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(splines[1]);
  ASSERT_THAT(spline1->GetKnotVector(1)->GetNumberOfKnots(), 6);
  ASSERT_THAT(spline2->GetKnotVector(1)->GetNumberOfKnots(), 7);
  ASSERT_THAT(spline1->GetKnotVector(1)->GetKnot(0).get(), DoubleEq(0));
  ASSERT_THAT(spline1->GetKnotVector(1)->GetKnot(2).get(), DoubleEq(0));
  ASSERT_THAT(spline1->GetKnotVector(1)->GetKnot(3).get(), DoubleEq(0.4));
  ASSERT_THAT(spline1->GetKnotVector(1)->GetKnot(5).get(), DoubleEq(0.4));
  ASSERT_THAT(spline2->GetKnotVector(1)->GetKnot(0).get(), DoubleEq(0.4));
  ASSERT_THAT(spline2->GetKnotVector(1)->GetKnot(2).get(), DoubleEq(0.4));
  ASSERT_THAT(spline2->GetKnotVector(1)->GetKnot(4).get(), DoubleEq(1));
  ASSERT_THAT(spline2->GetKnotVector(1)->GetKnot(6).get(), DoubleEq(1));
  for (int i = 0; i <= 100; ++i) {
    for (int j = 0; j <= 100; ++j) {
      std::array<ParamCoord, 2> param_coord{ParamCoord(i / 100.0), ParamCoord(j / 100.0)};
      if (j / 100.0 <= 0.4) {
        ASSERT_THAT(spline1->Evaluate(param_coord, {0})[0],
                    DoubleNear(bspline_2d_->Evaluate(param_coord, {0})[0], 0.000001));
      }
      if (j / 100.0 >= 0.4) {
        ASSERT_THAT(spline2->Evaluate(param_coord, {0})[0],
                    DoubleNear(bspline_2d_->Evaluate(param_coord, {0})[0], 0.000001));
      }
    }
  }
}

class NURBSFig5_9 : public Test {  // NOLINT
 public:
  NURBSFig5_9() {
    std::array<Degree, 2> degree = {Degree{3}, Degree{2}};
    std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector = {
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                           ParamCoord{1}, ParamCoord{1}, ParamCoord{1},
                                                           ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5},
                                                           ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({2.5, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.5, 0.0, 0.5})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 3.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 0.0, 3.0})),

        baf::ControlPoint(std::vector<double>({2.5, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({4.0, 3.0, 4.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 3.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 1.0, 3.0})),

        baf::ControlPoint(std::vector<double>({2.5, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.5, 0.5, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.5, 2.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 2.5, 2.0})),

        baf::ControlPoint(std::vector<double>({5.0, 5.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 4.5, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 4.0, 2.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 4.0, 2.0}))
    };
    std::vector<double> weights = {2, 4, 1, 1, 2, 0.5, 1, 1, 4, 4, 1, 1, 1, 1, 1, 1};
    nurbs_2d_ = std::make_shared<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<2>> nurbs_2d_;
};

TEST_F(NURBSFig5_9, IsSubdividedAtKnot0_4InSecondDirection) {  // NOLINT
  auto splines = nurbs_2d_->SudivideSpline(ParamCoord{0.4}, 1);
  auto spline1 = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(splines[0]);
  auto spline2 = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(splines[1]);
  ASSERT_THAT(spline1->GetKnotVector(1)->GetNumberOfKnots(), 6);
  ASSERT_THAT(spline2->GetKnotVector(1)->GetNumberOfKnots(), 7);
  ASSERT_THAT(spline1->GetKnotVector(1)->GetKnot(0).get(), DoubleEq(0));
  ASSERT_THAT(spline1->GetKnotVector(1)->GetKnot(2).get(), DoubleEq(0));
  ASSERT_THAT(spline1->GetKnotVector(1)->GetKnot(3).get(), DoubleEq(0.4));
  ASSERT_THAT(spline1->GetKnotVector(1)->GetKnot(5).get(), DoubleEq(0.4));
  ASSERT_THAT(spline2->GetKnotVector(1)->GetKnot(0).get(), DoubleEq(0.4));
  ASSERT_THAT(spline2->GetKnotVector(1)->GetKnot(2).get(), DoubleEq(0.4));
  ASSERT_THAT(spline2->GetKnotVector(1)->GetKnot(4).get(), DoubleEq(1));
  ASSERT_THAT(spline2->GetKnotVector(1)->GetKnot(6).get(), DoubleEq(1));
  for (int i = 0; i <= 100; ++i) {
    for (int j = 0; j <= 100; ++j) {
      std::array<ParamCoord, 2> param_coord{ParamCoord(i / 100.0), ParamCoord(j / 100.0)};
      if (j / 100.0 <= 0.4) {
        ASSERT_THAT(spline1->Evaluate(param_coord, {0})[0],
                    DoubleNear(nurbs_2d_->Evaluate(param_coord, {0})[0], 0.000001));
      }
      if (j / 100.0 >= 0.4) {
        ASSERT_THAT(spline2->Evaluate(param_coord, {0})[0],
                    DoubleNear(nurbs_2d_->Evaluate(param_coord, {0})[0], 0.000001));
      }
    }
  }
}

class Random3DBSpline : public Test {  // NOLINT
 public:
  Random3DBSpline() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0}, ParamCoord{1}};
    spl::RandomBSplineGenerator<3> spline_generator(limits, 10, 3);
    spl::BSpline<3> b_spline(*spline_generator.GetParameterSpace(), *spline_generator.GetPhysicalSpace());
    b_spline_3d_ = std::make_shared<spl::BSpline<3>>(b_spline);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> b_spline_3d_;
};

TEST_F(Random3DBSpline, IsSubdividedAtKnot0_4InSecondDirection) {  // NOLINT
  auto splines = b_spline_3d_->SudivideSpline(ParamCoord{0.4}, 1);
  auto spline1 = std::any_cast<std::shared_ptr<spl::BSpline<3>>>(splines[0]);
  auto spline2 = std::any_cast<std::shared_ptr<spl::BSpline<3>>>(splines[1]);
  double s = 15;
  for (int i = 0; i <= s; ++i) {
    for (int j = 0; j <= s; ++j) {
      for (int k = 0; k <= s; ++k) {
        std::array<ParamCoord, 3> param_coord{ParamCoord(i / s), ParamCoord(j / s), ParamCoord(k / s)};
        if (j / s <= 0.4) {
          ASSERT_THAT(spline1->Evaluate(param_coord, {0})[0],
                      DoubleNear(b_spline_3d_->Evaluate(param_coord, {0})[0], 0.000001));
        }
        if (j / s >= 0.4) {
          ASSERT_THAT(spline2->Evaluate(param_coord, {0})[0],
                      DoubleNear(b_spline_3d_->Evaluate(param_coord, {0})[0], 0.000001));
        }
      }
    }
  }
}
