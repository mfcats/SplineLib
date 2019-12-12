/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/spl/b_spline.h"
#include "src/spl/nurbs.h"
#include "test/spl/random/random_spline_utils.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

using namespace splinelib::src;

class BSplineFig5_9 : public Test {  // NOLINT
 public:
  BSplineFig5_9() {
    std::array<Degree, 2> degree = {Degree{3}, Degree{2}};
    baf::KnotVectors<2> knot_vector = {
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0},
             ParametricCoordinate{1},
             ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.5},
             ParametricCoordinate{1},
             ParametricCoordinate{1}, ParametricCoordinate{1}}))};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint(std::vector<double>({5.0, 0.0, 2.0})),
        spl::ControlPoint(std::vector<double>({2.0, 0.0, 2.0})),
        spl::ControlPoint(std::vector<double>({1.0, 0.0, 3.0})),
        spl::ControlPoint(std::vector<double>({-1.0, 0.0, 3.0})),

        spl::ControlPoint(std::vector<double>({5.0, 2.0, 2.0})),
        spl::ControlPoint(std::vector<double>({2.0, 1.5, 2.0})),
        spl::ControlPoint(std::vector<double>({1.0, 1.0, 3.0})),
        spl::ControlPoint(std::vector<double>({-1.0, 1.0, 3.0})),

        spl::ControlPoint(std::vector<double>({5.0, 2.0, 0.0})),
        spl::ControlPoint(std::vector<double>({2.0, 2.0, 0.0})),
        spl::ControlPoint(std::vector<double>({1.0, 2.5, 2.0})),
        spl::ControlPoint(std::vector<double>({-1.0, 2.5, 2.0})),

        spl::ControlPoint(std::vector<double>({5.0, 5.0, 0.0})),
        spl::ControlPoint(std::vector<double>({2.0, 4.5, 0.0})),
        spl::ControlPoint(std::vector<double>({1.0, 4.0, 2.0})),
        spl::ControlPoint(std::vector<double>({-1.0, 4.0, 2.0}))
    };
    bspline_2d_ = std::make_shared<spl::BSpline<2>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> bspline_2d_;
};

TEST_F(BSplineFig5_9, IsSubdividedAtKnot0_7InFirstDirection) {  // NOLINT
  auto splines = bspline_2d_->SudivideSpline(ParametricCoordinate{0.7}, 0);
  ASSERT_THAT(splines[0]->GetKnotVector(0)->GetNumberOfKnots(), 8);
  ASSERT_THAT(splines[1]->GetKnotVector(0)->GetNumberOfKnots(), 8);
  ASSERT_THAT(splines[0]->GetKnotVector(0)->GetFirstKnot().Get(), DoubleEq(0));
  ASSERT_THAT((*splines[0]->GetKnotVector(0))[3].Get(), DoubleEq(0));
  ASSERT_THAT((*splines[0]->GetKnotVector(0))[4].Get(), DoubleEq(0.7));
  ASSERT_THAT((*splines[0]->GetKnotVector(0))[7].Get(), DoubleEq(0.7));
  ASSERT_THAT(splines[1]->GetKnotVector(0)->GetFirstKnot().Get(), DoubleEq(0.7));
  ASSERT_THAT((*splines[1]->GetKnotVector(0))[3].Get(), DoubleEq(0.7));
  ASSERT_THAT((*splines[1]->GetKnotVector(0))[4].Get(), DoubleEq(1));
  ASSERT_THAT((*splines[1]->GetKnotVector(0))[7].Get(), DoubleEq(1));
  double steps = 100;
  for (int i = 0; i <= steps; ++i) {
    ParametricCoordinate coord2 = ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)};
    std::array<ParametricCoordinate, 2> param_coord{ParametricCoordinate(i / steps), coord2};
    int spline_number = i / steps >= 0.7 ? 1 : 0;
    std::vector<double> evaluated_splitted_spline = splines[spline_number]->Evaluate(param_coord, {0, 1});
    std::vector<double> original_spline = bspline_2d_->Evaluate(param_coord, {0, 1});
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(evaluated_splitted_spline[j], DoubleNear(original_spline[j], 0.000001));
    }
  }
}

class Random1DBSplineToSplit : public Test {  // NOLINT
 public:
  Random1DBSplineToSplit() {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0}, ParametricCoordinate{1}};
    b_spline_1d_ = splinelib::test::RandomSplineUtils<1>::GenerateRandomBSpline(limits, 10, 2);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> b_spline_1d_;
};

TEST_F(Random1DBSplineToSplit, IsSubdividedAtKnot0_25) {  // NOLINT
  auto splines = b_spline_1d_->SudivideSpline(ParametricCoordinate{0.25}, 0);
  double steps = 100;
  for (int i = 0; i <= steps; ++i) {
    std::array<ParametricCoordinate, 1> param_coord{ParametricCoordinate(i / steps)};
    int spline_number = i / steps >= 0.25 ? 1 : 0;
    double evaluated_splitted_spline = splines[spline_number]->Evaluate(param_coord, {0})[0];
    double original_spline = b_spline_1d_->Evaluate(param_coord, {0})[0];
    ASSERT_THAT(evaluated_splitted_spline, DoubleNear(original_spline, 0.000001));
  }
}

class Random1DNURBSToSplit : public Test {  // NOLINT
 public:
  Random1DNURBSToSplit() {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0}, ParametricCoordinate{1}};
    nurbs_1d_ = splinelib::test::RandomSplineUtils<1>::GenerateRandomNURBS(limits, 10, 3);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_;
};

TEST_F(Random1DNURBSToSplit, IsSubdividedAtKnot0_99) {  // NOLINT
  auto splines = nurbs_1d_->SudivideSpline(ParametricCoordinate{0.99}, 0);
  double steps = 100;
  for (int i = 0; i <= steps; ++i) {
    std::array<ParametricCoordinate, 1> param_coord{ParametricCoordinate(i / steps)};
    int spline_number = i / steps >= 0.99 ? 1 : 0;
    double evaluated_splitted_spline = splines[spline_number]->Evaluate(param_coord, {0})[0];
    double original_spline = nurbs_1d_->Evaluate(param_coord, {0})[0];
    ASSERT_THAT(evaluated_splitted_spline, DoubleNear(original_spline, 0.000001));
  }
}

class Random2DBSplineToSplit : public Test {  // NOLINT
 public:
  Random2DBSplineToSplit() {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0}, ParametricCoordinate{1}};
    b_spline_2d_ = splinelib::test::RandomSplineUtils<2>::GenerateRandomBSpline(limits, 10, 3);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> b_spline_2d_;
};

TEST_F(Random2DBSplineToSplit, IsSubdividedAtKnot0_1InFirstDirection) {  // NOLINT
  auto splines = b_spline_2d_->SudivideSpline(ParametricCoordinate{0.1}, 0);
  double steps = 50;
  for (int i = 0; i <= steps; ++i) {
    ParametricCoordinate coord2 = ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)};
    std::array<ParametricCoordinate, 2> param_coord{ParametricCoordinate(i / steps), coord2};
    int spline_number = i / steps >= 0.1 ? 1 : 0;
    std::vector<double> evaluated_splitted_spline = splines[spline_number]->Evaluate(param_coord, {0, 1});
    std::vector<double> original_spline = b_spline_2d_->Evaluate(param_coord, {0, 1});
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(evaluated_splitted_spline[j], DoubleNear(original_spline[j], 0.000001));
    }
  }
}

class Random2DNURBSToSplit : public Test {  // NOLINT
 public:
  Random2DNURBSToSplit() {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0}, ParametricCoordinate{1}};
    nurbs_2d_ = splinelib::test::RandomSplineUtils<2>::GenerateRandomNURBS(limits, 10, 3);
  }

 protected:
  std::shared_ptr<spl::NURBS<2>> nurbs_2d_;
};

TEST_F(Random2DNURBSToSplit, IsSubdividedAtKnot0_5InSecondDirection) {  // NOLINT
  auto splines = nurbs_2d_->SudivideSpline(ParametricCoordinate{0.5}, 1);
  double steps = 50;
  for (int i = 0; i <= steps; ++i) {
    ParametricCoordinate coord1 = ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)};
    std::array<ParametricCoordinate, 2> param_coord{coord1, ParametricCoordinate(i / steps)};
    int spline_number = i / steps >= 0.5 ? 1 : 0;
    std::vector<double> evaluated_splitted_spline = splines[spline_number]->Evaluate(param_coord, {0, 1});
    std::vector<double> original_spline = nurbs_2d_->Evaluate(param_coord, {0, 1});
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(evaluated_splitted_spline[j], DoubleNear(original_spline[j], 0.000001));
    }
  }
}

class Random3DBSplineToSplit : public Test {  // NOLINT
 public:
  Random3DBSplineToSplit() {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0}, ParametricCoordinate{1}};
    b_spline_3d_ = splinelib::test::RandomSplineUtils<3>::GenerateRandomBSpline(limits, 10, 3);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> b_spline_3d_;
};

TEST_F(Random3DBSplineToSplit, IsSubdividedAtKnot0_4InSecondDirection) {  // NOLINT
  auto splines = b_spline_3d_->SudivideSpline(ParametricCoordinate{0.4}, 1);
  double steps = 25;
  for (int i = 0; i <= steps; ++i) {
    ParametricCoordinate coord1 = ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)};
    ParametricCoordinate coord3 = ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)};
    std::array<ParametricCoordinate, 3> param_coord{coord1, ParametricCoordinate(i / steps), coord3};
    int spline_number = i / steps >= 0.4 ? 1 : 0;
    std::vector<double> evaluated_splitted_spline = splines[spline_number]->Evaluate(param_coord, {0, 1, 2});
    std::vector<double> original_spline = b_spline_3d_->Evaluate(param_coord, {0, 1, 2});
    for (int j = 0; j < 3; ++j) {
      ASSERT_THAT(evaluated_splitted_spline[j], DoubleNear(original_spline[j], 0.000001));
    }
  }
}
