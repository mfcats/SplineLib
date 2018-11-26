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

class BSplineEx5_1 : public Test {  // NOLINT
 public:
  BSplineEx5_1() {
    std::array<Degree, 1> degree = {Degree{3}};
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2},
                         ParamCoord{3}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}}))};
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_after = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2},
                         ParamCoord{3}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({4.0, 0.0})),
        baf::ControlPoint(std::vector<double>({5.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 4.0})),
        baf::ControlPoint(std::vector<double>({2.0, 4.0})),
        baf::ControlPoint(std::vector<double>({1.3, 2.3}))
    };
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(knot_vector_after, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(BSplineEx5_1, InsertsKnot2_5Correctly) {  // NOLINT
  bspline_1d_after_->InsertKnot(ParamCoord(2.5), 0);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(6).get(), DoubleEq(2.5));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() + 1);
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint(std::vector<double>({0.0, 1.0})),
      baf::ControlPoint(std::vector<double>({1.0, 2.0})),
      baf::ControlPoint(std::vector<double>({2.0, 0.0})),
      baf::ControlPoint(std::vector<double>({11.0 / 3.0, 0.0})),
      baf::ControlPoint(std::vector<double>({4.5, 1.0})),
      baf::ControlPoint(std::vector<double>({29.0 / 6.0, 7.0 / 3.0})),
      baf::ControlPoint(std::vector<double>({4.0, 4.0})),
      baf::ControlPoint(std::vector<double>({2.0, 4.0})),
      baf::ControlPoint(std::vector<double>({1.3, 2.3}))
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, j), DoubleEq(new_control_points[i].GetValue(j)));
    }
  }
  for (int i = 0; i <= 50; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / 10.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleEq(bspline_1d_before_->Evaluate(param_coord, {0})[0]));
  }
}

class NURBSEx5_2 : public Test {  // NOLINT
 public:
  NURBSEx5_2() {
    std::array<Degree, 1> degree = {Degree{3}};
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2},
                         ParamCoord{3}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}}))};
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_after = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2},
                         ParamCoord{3}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({5.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.3, 2.3}))
    };
    std::vector<double> weights = {1.0, 1.0, 2.0, 4.0, 1.0, 4.0, 2.0, 1.0};
    nurbs_1d_before_ = std::make_shared<spl::NURBS<1>>(knot_vector_before, degree, control_points, weights);
    nurbs_1d_after_ = std::make_shared<spl::NURBS<1>>(knot_vector_after, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_before_;
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_after_;
};

TEST_F(NURBSEx5_2, InsertsKnot2_0Correctly) {  // NOLINT
  nurbs_1d_after_->InsertKnot(ParamCoord(2.0), 0);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetKnot(6).get(), DoubleEq(2.0));
  ASSERT_THAT(nurbs_1d_after_->GetNumberOfControlPoints(), nurbs_1d_before_->GetNumberOfControlPoints() + 1);
  ASSERT_THAT(nurbs_1d_after_->GetWeights().size(), nurbs_1d_after_->GetNumberOfControlPoints());
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint(std::vector<double>({0.0, 1.0})),
      baf::ControlPoint(std::vector<double>({1.0, 2.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({13.0 / 9.0, 2.0 / 9.0})),
      baf::ControlPoint(std::vector<double>({5.0, 2.0})),
      baf::ControlPoint(std::vector<double>({1.0, 1.0})),
      baf::ControlPoint(std::vector<double>({1.0, 2.0})),
      baf::ControlPoint(std::vector<double>({1.3, 2.3}))
  };
  for (int i = 3; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(nurbs_1d_after_->GetControlPoint({i}, j), DoubleEq(new_control_points[i].GetValue(j)));
    }
  }
  ASSERT_THAT(nurbs_1d_after_->GetWeights().size(), nurbs_1d_before_->GetWeights().size() + 1);
  std::vector<double> new_weights = {1.0, 1.0, 2.0, 10.0 / 3.0, 3.0, 1.0, 4.0, 2.0, 1.0};
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    ASSERT_THAT(nurbs_1d_after_->GetWeight({i}), DoubleEq(new_weights[i]));
  }
  for (int i = 0; i <= 50; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / 10.0)};
    ASSERT_THAT(nurbs_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleEq(nurbs_1d_before_->Evaluate(param_coord, {0})[0]));
  }
}

class Random2DBSplineForKnotInsertion : public Test {  // NOLINT
 public:
  Random2DBSplineForKnotInsertion() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomBSplineGenerator<2> spline_generator(limits, 10, 3);
    spl::BSpline<2> b_spline(*spline_generator.GetParameterSpace(), *spline_generator.GetPhysicalSpace());
    bspline_2d_before_ = std::make_shared<spl::BSpline<2>>(b_spline);
    spl::BSpline<2> b_spline_after(b_spline);
    bspline_2d_after_ = std::make_shared<spl::BSpline<2>>(b_spline_after);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> bspline_2d_before_;
  std::shared_ptr<spl::BSpline<2>> bspline_2d_after_;
};

TEST_F(Random2DBSplineForKnotInsertion, InsertsKnot0_4InFirstDirectionCorrectly) {  // NOLINT
  bspline_2d_after_->InsertKnot(ParamCoord(0.4), 0);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(1)->GetNumberOfKnots());
  ASSERT_THAT(bspline_2d_after_->GetPointsPerDirection()[0], bspline_2d_before_->GetPointsPerDirection()[0] + 1);
  double steps = 50;
  for (int i = 0; i <= steps; ++i) {
    ParamCoord coord2 = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    std::array<ParamCoord, 2> param_coord{ParamCoord(i / steps), coord2};
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(bspline_2d_after_->Evaluate(param_coord, {j})[0],
                  DoubleNear(bspline_2d_before_->Evaluate(param_coord, {j})[0], 0.000001));
    }
  }
}

TEST_F(Random2DBSplineForKnotInsertion, InsertsKnot0_7InSecondDirectionCorrectly) {  // NOLINT
  bspline_2d_after_->InsertKnot(ParamCoord(0.7), 1);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(1)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_2d_after_->GetPointsPerDirection()[1], bspline_2d_before_->GetPointsPerDirection()[1] + 1);
  double steps = 50;
  for (int i = 0; i <= steps; ++i) {
    ParamCoord coord1 = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    std::array<ParamCoord, 2> param_coord{coord1, ParamCoord(i / steps)};
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(bspline_2d_after_->Evaluate(param_coord, {j})[0],
                  DoubleNear(bspline_2d_before_->Evaluate(param_coord, {j})[0], 0.000001));
    }
  }
}

TEST_F(Random2DBSplineForKnotInsertion, InsertsKnot0_4InFirstAndKnot0_7InSecondDirectionCorrectly) {  // NOLINT
  bspline_2d_after_->InsertKnot(ParamCoord(0.4), 0);
  bspline_2d_after_->InsertKnot(ParamCoord(0.7), 1);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(1)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_2d_after_->GetPointsPerDirection()[0], bspline_2d_before_->GetPointsPerDirection()[0] + 1);
  ASSERT_THAT(bspline_2d_after_->GetPointsPerDirection()[1], bspline_2d_before_->GetPointsPerDirection()[1] + 1);
  double steps = 10;
  for (int i = 0; i <= steps; ++i) {
    for (int j = 0; j <= steps; ++j) {
      std::array<ParamCoord, 2> param_coord{ParamCoord(i / steps), ParamCoord(j / steps)};
      for (int k = 0; k < 2; ++k) {
        ASSERT_THAT(bspline_2d_after_->Evaluate(param_coord, {k})[0],
                    DoubleNear(bspline_2d_before_->Evaluate(param_coord, {k})[0], 0.000001));
      }
    }
  }
}

class BSpline3DEx : public Test {  // NOLINT
 public:
  BSpline3DEx() {
    std::array<Degree, 3> degree = {Degree{2}, Degree{1}, Degree{2}};
    std::array<std::shared_ptr<baf::KnotVector>, 3> knot_vector_before = {
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}, ParamCoord{1}}))};
    std::array<std::shared_ptr<baf::KnotVector>, 3> knot_vector_after = {
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0, 2.0}))
    };
    bspline_3d_before_ = std::make_shared<spl::BSpline<3>>(knot_vector_before, degree, control_points);
    bspline_3d_after_ = std::make_shared<spl::BSpline<3>>(knot_vector_after, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> bspline_3d_before_;
  std::shared_ptr<spl::BSpline<3>> bspline_3d_after_;
};

TEST_F(BSpline3DEx, InsertsKnot0_4InFirst_Knot0_99InSecond_Knot0_01InThirdDirectionCorrectly) {  // NOLINT
  bspline_3d_after_->InsertKnot(ParamCoord(0.4), 0);
  bspline_3d_after_->InsertKnot(ParamCoord(0.99), 1);
  bspline_3d_after_->InsertKnot(ParamCoord(0.01), 2);
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_3d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              bspline_3d_before_->GetKnotVector(1)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(2)->GetNumberOfKnots(),
              bspline_3d_before_->GetKnotVector(2)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(0)->GetKnot(3).get(), DoubleEq(0.4));
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(1)->GetKnot(2).get(), DoubleEq(0.99));
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(2)->GetKnot(3).get(), DoubleEq(0.01));
  ASSERT_THAT(bspline_3d_after_->GetNumberOfControlPoints(), bspline_3d_before_->GetNumberOfControlPoints() + 30);
  ASSERT_THAT(bspline_3d_after_->GetPointsPerDirection()[0], bspline_3d_before_->GetPointsPerDirection()[0] + 1);
  ASSERT_THAT(bspline_3d_after_->GetPointsPerDirection()[1], bspline_3d_before_->GetPointsPerDirection()[1] + 1);
  ASSERT_THAT(bspline_3d_after_->GetPointsPerDirection()[2], bspline_3d_before_->GetPointsPerDirection()[2] + 1);
  double s = 20;
  for (int i = 0; i <= s; ++i) {
    for (int j = 0; j <= s; ++j) {
      for (int k = 0; k <= s; ++k) {
        std::array<ParamCoord, 3> param_coord{ParamCoord(i / s), ParamCoord(j / s), ParamCoord(k / s)};
        ASSERT_THAT(bspline_3d_after_->Evaluate(param_coord, {0})[0],
                    DoubleNear(bspline_3d_before_->Evaluate(param_coord, {0})[0], 0.000001));
      }
    }
  }
}

class NURBS3DEx : public Test {  // NOLINT
 public:
  NURBS3DEx() {
    std::array<Degree, 3> degree = {Degree{2}, Degree{1}, Degree{2}};
    std::array<std::shared_ptr<baf::KnotVector>, 3> knot_vector_before = {
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}, ParamCoord{1}}))};
    std::array<std::shared_ptr<baf::KnotVector>, 3> knot_vector_after = {
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                           ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0, 2.0}))
    };
    std::vector<double>
        weights = {0.3, 0.5, 0.8, 1.0, 0.9, 1.1, 2.0, 1.5, 1.4, 1.9, 2.2, 2.1, 3.0, 3.1, 4.0, 2.9, 3.8, 4.2};
    nurbs_3d_before_ = std::make_shared<spl::NURBS<3>>(knot_vector_before, degree, control_points, weights);
    nurbs_3d_after_ = std::make_shared<spl::NURBS<3>>(knot_vector_after, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<3>> nurbs_3d_before_;
  std::shared_ptr<spl::NURBS<3>> nurbs_3d_after_;
};

TEST_F(NURBS3DEx, InsertsKnot0_4InFirst_Knot0_99InSecond_Knot0_01InThirdDirectionCorrectly) {  // NOLINT
  nurbs_3d_after_->InsertKnot(ParamCoord(0.4), 0);
  nurbs_3d_after_->InsertKnot(ParamCoord(0.99), 1);
  nurbs_3d_after_->InsertKnot(ParamCoord(0.01), 2);
  ASSERT_THAT(nurbs_3d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_3d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(nurbs_3d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              nurbs_3d_before_->GetKnotVector(1)->GetNumberOfKnots() + 1);
  ASSERT_THAT(nurbs_3d_after_->GetKnotVector(2)->GetNumberOfKnots(),
              nurbs_3d_before_->GetKnotVector(2)->GetNumberOfKnots() + 1);
  ASSERT_THAT(nurbs_3d_after_->GetKnotVector(0)->GetKnot(3).get(), DoubleEq(0.4));
  ASSERT_THAT(nurbs_3d_after_->GetKnotVector(1)->GetKnot(2).get(), DoubleEq(0.99));
  ASSERT_THAT(nurbs_3d_after_->GetKnotVector(2)->GetKnot(3).get(), DoubleEq(0.01));
  ASSERT_THAT(nurbs_3d_after_->GetNumberOfControlPoints(), nurbs_3d_before_->GetNumberOfControlPoints() + 30);
  ASSERT_THAT(nurbs_3d_after_->GetWeights().size(), nurbs_3d_after_->GetNumberOfControlPoints());
  ASSERT_THAT(nurbs_3d_after_->GetPointsPerDirection()[0], nurbs_3d_before_->GetPointsPerDirection()[0] + 1);
  ASSERT_THAT(nurbs_3d_after_->GetPointsPerDirection()[1], nurbs_3d_before_->GetPointsPerDirection()[1] + 1);
  ASSERT_THAT(nurbs_3d_after_->GetPointsPerDirection()[2], nurbs_3d_before_->GetPointsPerDirection()[2] + 1);
  double s = 15;
  for (int i = 0; i <= s; ++i) {
    for (int j = 0; j <= s; ++j) {
      for (int k = 0; k <= s; ++k) {
        std::array<ParamCoord, 3> param_coord{ParamCoord(i / s), ParamCoord(j / s), ParamCoord(k / s)};
        ASSERT_THAT(nurbs_3d_after_->Evaluate(param_coord, {0})[0],
                    DoubleNear(nurbs_3d_before_->Evaluate(param_coord, {0})[0], 0.000001));
      }
    }
  }
}

class BSplineFig5_16 : public Test {  // NOLINT
 public:
  BSplineFig5_16() {
    std::array<Degree, 1> degree = {Degree{3}};
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.7},
                         ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_after = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.7},
                         ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.15})),
        baf::ControlPoint(std::vector<double>({0.5, 1.5})),
        baf::ControlPoint(std::vector<double>({2.3, 2.0})),
        baf::ControlPoint(std::vector<double>({1.7, 0.3})),
        baf::ControlPoint(std::vector<double>({3.0, 0.0}))
    };
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(knot_vector_after, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(BSplineFig5_16, InsertsMidpoints) {  // NOLINT
  std::vector<ParamCoord> new_knots = {ParamCoord{0.15}, ParamCoord{0.5}, ParamCoord{0.85}};
  bspline_1d_after_->RefineKnots(new_knots, 0);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 3);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(0.15));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(5).get(), DoubleEq(0.3));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(6).get(), DoubleEq(0.5));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(7).get(), DoubleEq(0.7));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(8).get(), DoubleEq(0.85));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() + 3);
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint(std::vector<double>({0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.5, 0.075})),
      baf::ControlPoint(std::vector<double>({0.892857, 0.439286})),
      baf::ControlPoint(std::vector<double>({0.805102, 1.25051})),
      baf::ControlPoint(std::vector<double>({1.4, 1.75})),
      baf::ControlPoint(std::vector<double>({1.97245, 1.5648})),
      baf::ControlPoint(std::vector<double>({1.82857, 0.664286})),
      baf::ControlPoint(std::vector<double>({2.35, 0.15})),
      baf::ControlPoint(std::vector<double>({3.0, 0.0}))
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(bspline_1d_after_->GetControlPoint({i}, j), DoubleNear(new_control_points[i].GetValue(j), 0.00001));
    }
  }
  for (int i = 0; i <= 100; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / 100.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleNear(bspline_1d_before_->Evaluate(param_coord, {0})[0], 0.00001));
  }
}

TEST_F(BSplineFig5_16, InsertsKnot0_5MultipleTimesWithKnotRefinement) {  // NOLINT
  std::vector<ParamCoord> new_knots = {ParamCoord{0.5}, ParamCoord{0.5}, ParamCoord{0.5}};
  bspline_1d_after_->RefineKnots(new_knots, 0);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 3);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(0.3));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(5).get(), DoubleEq(0.5));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(6).get(), DoubleEq(0.5));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(7).get(), DoubleEq(0.5));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(8).get(), DoubleEq(0.7));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() + 3);
  for (int i = 0; i <= 100; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / 100.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleNear(bspline_1d_before_->Evaluate(param_coord, {0})[0], 0.00001));
  }
}

TEST_F(BSplineFig5_16, InsertsKnot0_5MultipleTimes) {  // NOLINT
  bspline_1d_after_->InsertKnot(ParamCoord{0.5}, 0, 3);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 3);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(0.3));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(5).get(), DoubleEq(0.5));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(6).get(), DoubleEq(0.5));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(7).get(), DoubleEq(0.5));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(8).get(), DoubleEq(0.7));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() + 3);
  for (int i = 0; i <= 100; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / 100.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleNear(bspline_1d_before_->Evaluate(param_coord, {0})[0], 0.00001));
  }
}

TEST_F(BSplineFig5_16, InsertsKnot0_3MultipleTimes) {  // NOLINT
  bspline_1d_after_->InsertKnot(ParamCoord{0.3}, 0, 2);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 2);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(0.3));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(5).get(), DoubleEq(0.3));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(6).get(), DoubleEq(0.3));
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(7).get(), DoubleEq(0.7));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() + 2);
  for (int i = 0; i <= 100; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / 100.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleNear(bspline_1d_before_->Evaluate(param_coord, {0})[0], 0.00001));
  }
}
