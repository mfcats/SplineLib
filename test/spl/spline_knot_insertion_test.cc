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

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

using namespace splinelib::src;

class BSpline1DEx5_1 : public Test {  // NOLINT
 public:
  BSpline1DEx5_1() {
    std::array<Degree, 1> degree = {Degree{3}};
    baf::KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(baf::KnotVector(
        {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0},
         ParametricCoordinate{1},
         ParametricCoordinate{2}, ParametricCoordinate{3}, ParametricCoordinate{4}, ParametricCoordinate{5},
         ParametricCoordinate{5},
         ParametricCoordinate{5}, ParametricCoordinate{5}}))};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint(std::vector<double>({0.0, 1.0})),
        spl::ControlPoint(std::vector<double>({1.0, 2.0})),
        spl::ControlPoint(std::vector<double>({2.0, 0.0})),
        spl::ControlPoint(std::vector<double>({4.0, 0.0})),
        spl::ControlPoint(std::vector<double>({5.0, 2.0})),
        spl::ControlPoint(std::vector<double>({4.0, 4.0})),
        spl::ControlPoint(std::vector<double>({2.0, 4.0})),
        spl::ControlPoint(std::vector<double>({1.3, 2.3}))
    };
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    spl::BSpline<1> b_spline_after(*bspline_1d_before_);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(b_spline_after);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(BSpline1DEx5_1, InsertsKnot2_5Correctly) {  // NOLINT
  bspline_1d_after_->InsertKnot(ParametricCoordinate(2.5), 0);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[6].Get(), DoubleEq(2.5));
  ASSERT_THAT(bspline_1d_after_->GetTotalNumberOfControlPoints(),
              bspline_1d_before_->GetTotalNumberOfControlPoints() + 1);
  std::vector<spl::ControlPoint> new_control_points = {
      spl::ControlPoint(std::vector<double>({0.0, 1.0})),
      spl::ControlPoint(std::vector<double>({1.0, 2.0})),
      spl::ControlPoint(std::vector<double>({2.0, 0.0})),
      spl::ControlPoint(std::vector<double>({11.0 / 3.0, 0.0})),
      spl::ControlPoint(std::vector<double>({4.5, 1.0})),
      spl::ControlPoint(std::vector<double>({29.0 / 6.0, 7.0 / 3.0})),
      spl::ControlPoint(std::vector<double>({4.0, 4.0})),
      spl::ControlPoint(std::vector<double>({2.0, 4.0})),
      spl::ControlPoint(std::vector<double>({1.3, 2.3}))
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(bspline_1d_after_->GetControlPoint(i, j),
          DoubleEq(new_control_points[i][Dimension{j}]));
    }
  }
  for (int i = 0; i <= 50; ++i) {
    std::array<ParametricCoordinate, 1> param_coord{ParametricCoordinate(i / 10.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleEq(bspline_1d_before_->Evaluate(param_coord, {0})[0]));
  }
}

class NURBS1DEx5_2 : public Test {  // NOLINT
 public:
  NURBS1DEx5_2() {
    std::array<Degree, 1> degree = {Degree{3}};
    baf::KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(baf::KnotVector(
        {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0},
         ParametricCoordinate{1},
         ParametricCoordinate{2}, ParametricCoordinate{3}, ParametricCoordinate{4}, ParametricCoordinate{5},
         ParametricCoordinate{5},
         ParametricCoordinate{5}, ParametricCoordinate{5}}))};
    baf::KnotVectors<1> knot_vector_after = {std::make_shared<baf::KnotVector>(baf::KnotVector(
        {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0},
         ParametricCoordinate{1},
         ParametricCoordinate{2}, ParametricCoordinate{3}, ParametricCoordinate{4}, ParametricCoordinate{5},
         ParametricCoordinate{5},
         ParametricCoordinate{5}, ParametricCoordinate{5}}))};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint(std::vector<double>({0.0, 1.0})),
        spl::ControlPoint(std::vector<double>({1.0, 2.0})),
        spl::ControlPoint(std::vector<double>({1.0, 0.0})),
        spl::ControlPoint(std::vector<double>({1.0, 0.0})),
        spl::ControlPoint(std::vector<double>({5.0, 2.0})),
        spl::ControlPoint(std::vector<double>({1.0, 1.0})),
        spl::ControlPoint(std::vector<double>({1.0, 2.0})),
        spl::ControlPoint(std::vector<double>({1.3, 2.3}))
    };
    std::vector<Weight> weights = {Weight{1.0}, Weight{1.0}, Weight{2.0}, Weight{4.0}, Weight{1.0}, Weight{4.0},
                                   Weight{2.0}, Weight{1.0}};
    nurbs_1d_before_ = std::make_shared<spl::NURBS<1>>(knot_vector_before, degree, control_points, weights);
    spl::NURBS<1> nurbs_after(*nurbs_1d_before_);
    nurbs_1d_after_ = std::make_shared<spl::NURBS<1>>(nurbs_after);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_before_;
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_after_;
};

TEST_F(NURBS1DEx5_2, InsertsKnot2_0Correctly) {  // NOLINT
  nurbs_1d_after_->InsertKnot(ParametricCoordinate(2.0), 0);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT((*nurbs_1d_after_->GetKnotVector(0))[6].Get(), DoubleEq(2.0));
  ASSERT_THAT(nurbs_1d_after_->GetTotalNumberOfControlPoints(), nurbs_1d_before_->GetTotalNumberOfControlPoints() + 1);
  std::vector<spl::ControlPoint> new_control_points = {
      spl::ControlPoint(std::vector<double>({0.0, 1.0})),
      spl::ControlPoint(std::vector<double>({1.0, 2.0})),
      spl::ControlPoint(std::vector<double>({1.0, 0.0})),
      spl::ControlPoint(std::vector<double>({1.0, 0.0})),
      spl::ControlPoint(std::vector<double>({13.0 / 9.0, 2.0 / 9.0})),
      spl::ControlPoint(std::vector<double>({5.0, 2.0})),
      spl::ControlPoint(std::vector<double>({1.0, 1.0})),
      spl::ControlPoint(std::vector<double>({1.0, 2.0})),
      spl::ControlPoint(std::vector<double>({1.3, 2.3}))
  };
  for (int i = 3; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(nurbs_1d_after_->GetControlPoint(i, j),
          DoubleEq(new_control_points[i][Dimension{j}]));
    }
  }
  std::vector<double> new_weights = {1.0, 1.0, 2.0, 10.0 / 3.0, 3.0, 1.0, 4.0, 2.0, 1.0};
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    ASSERT_THAT(nurbs_1d_after_->GetWeight(i), DoubleEq(new_weights[i]));
  }
  for (int i = 0; i <= 50; ++i) {
    std::array<ParametricCoordinate, 1> param_coord{ParametricCoordinate(i / 10.0)};
    ASSERT_THAT(nurbs_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleEq(nurbs_1d_before_->Evaluate(param_coord, {0})[0]));
  }
}

class BSpline1DFig5_16 : public Test {  // NOLINT
 public:
  BSpline1DFig5_16() {
    std::array<Degree, 1> degree = {Degree{3}};
    baf::KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(baf::KnotVector(
        {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0},
         ParametricCoordinate{0.3},
         ParametricCoordinate{0.7}, ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{1},
         ParametricCoordinate{1}}))};
    baf::KnotVectors<1> knot_vector_after = {std::make_shared<baf::KnotVector>(baf::KnotVector(
        {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0},
         ParametricCoordinate{0.3},
         ParametricCoordinate{0.7}, ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{1},
         ParametricCoordinate{1}}))};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint(std::vector<double>({0.0, 0.0})),
        spl::ControlPoint(std::vector<double>({1.0, 0.15})),
        spl::ControlPoint(std::vector<double>({0.5, 1.5})),
        spl::ControlPoint(std::vector<double>({2.3, 2.0})),
        spl::ControlPoint(std::vector<double>({1.7, 0.3})),
        spl::ControlPoint(std::vector<double>({3.0, 0.0}))
    };
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(knot_vector_after, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(BSpline1DFig5_16, InsertsMidpoints) {  // NOLINT
  std::vector<ParametricCoordinate>
      new_knots = {ParametricCoordinate{0.15}, ParametricCoordinate{0.5}, ParametricCoordinate{0.85}};
  bspline_1d_after_->RefineKnots(new_knots, 0);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 3);
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[4].Get(), DoubleEq(0.15));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[5].Get(), DoubleEq(0.3));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[6].Get(), DoubleEq(0.5));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[7].Get(), DoubleEq(0.7));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[8].Get(), DoubleEq(0.85));
  ASSERT_THAT(bspline_1d_after_->GetTotalNumberOfControlPoints(),
              bspline_1d_before_->GetTotalNumberOfControlPoints() + 3);
  std::vector<spl::ControlPoint> new_control_points = {
      spl::ControlPoint(std::vector<double>({0.0, 0.0})),
      spl::ControlPoint(std::vector<double>({0.5, 0.075})),
      spl::ControlPoint(std::vector<double>({0.892857, 0.439286})),
      spl::ControlPoint(std::vector<double>({0.805102, 1.25051})),
      spl::ControlPoint(std::vector<double>({1.4, 1.75})),
      spl::ControlPoint(std::vector<double>({1.97245, 1.5648})),
      spl::ControlPoint(std::vector<double>({1.82857, 0.664286})),
      spl::ControlPoint(std::vector<double>({2.35, 0.15})),
      spl::ControlPoint(std::vector<double>({3.0, 0.0}))
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(bspline_1d_after_->GetControlPoint(i, j),
          DoubleNear(new_control_points[i][Dimension{j}], 0.00001));
    }
  }
  for (int i = 0; i <= 100; ++i) {
    std::array<ParametricCoordinate, 1> param_coord{ParametricCoordinate(i / 100.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleNear(bspline_1d_before_->Evaluate(param_coord, {0})[0], 0.00001));
  }
}

TEST_F(BSpline1DFig5_16, InsertsKnot0_5MultipleTimesWithKnotRefinement) {  // NOLINT
  std::vector<ParametricCoordinate>
      new_knots = {ParametricCoordinate{0.5}, ParametricCoordinate{0.5}, ParametricCoordinate{0.5}};
  bspline_1d_after_->RefineKnots(new_knots, 0);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 3);
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[4].Get(), DoubleEq(0.3));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[5].Get(), DoubleEq(0.5));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[6].Get(), DoubleEq(0.5));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[7].Get(), DoubleEq(0.5));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[8].Get(), DoubleEq(0.7));
  ASSERT_THAT(bspline_1d_after_->GetTotalNumberOfControlPoints(),
              bspline_1d_before_->GetTotalNumberOfControlPoints() + 3);
  for (int i = 0; i <= 100; ++i) {
    std::array<ParametricCoordinate, 1> param_coord{ParametricCoordinate(i / 100.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleNear(bspline_1d_before_->Evaluate(param_coord, {0})[0], 0.00001));
  }
}

TEST_F(BSpline1DFig5_16, InsertsKnot0_5MultipleTimes) {  // NOLINT
  bspline_1d_after_->InsertKnot(ParametricCoordinate{0.5}, 0, 3);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 3);
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[4].Get(), DoubleEq(0.3));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[5].Get(), DoubleEq(0.5));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[6].Get(), DoubleEq(0.5));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[7].Get(), DoubleEq(0.5));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[8].Get(), DoubleEq(0.7));
  ASSERT_THAT(bspline_1d_after_->GetTotalNumberOfControlPoints(),
              bspline_1d_before_->GetTotalNumberOfControlPoints() + 3);
  for (int i = 0; i <= 100; ++i) {
    std::array<ParametricCoordinate, 1> param_coord{ParametricCoordinate(i / 100.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleNear(bspline_1d_before_->Evaluate(param_coord, {0})[0], 0.00001));
  }
}

TEST_F(BSpline1DFig5_16, InsertsKnot0_3MultipleTimes) {  // NOLINT
  bspline_1d_after_->InsertKnot(ParametricCoordinate{0.3}, 0, 2);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 2);
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[4].Get(), DoubleEq(0.3));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[5].Get(), DoubleEq(0.3));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[6].Get(), DoubleEq(0.3));
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[7].Get(), DoubleEq(0.7));
  ASSERT_THAT(bspline_1d_after_->GetTotalNumberOfControlPoints(),
              bspline_1d_before_->GetTotalNumberOfControlPoints() + 2);
  for (int i = 0; i <= 100; ++i) {
    std::array<ParametricCoordinate, 1> param_coord{ParametricCoordinate(i / 100.0)};
    ASSERT_THAT(bspline_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleNear(bspline_1d_before_->Evaluate(param_coord, {0})[0], 0.00001));
  }
}
