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

using namespace splinelib::src;

class BSpline1DFig5_26 : public Test {  // NOLINT
 public:
  BSpline1DFig5_26() {
    std::array<Degree, 1> degree = {Degree{3}};
    ParametricCoordinate zero(0), one(1), two(2);
    baf::KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({zero, zero, zero, zero, one, one, one, two, two, two, two}))};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint({0.0, 0.0}), spl::ControlPoint({0.0, 1.5}), spl::ControlPoint({1.0, 2.0}),
        spl::ControlPoint({2.0, 2.0}), spl::ControlPoint({3.0, 2.0}), spl::ControlPoint({4.0, 1.5}),
        spl::ControlPoint({4.0, 0.0})
    };
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    spl::BSpline<1> b_spline_after(*bspline_1d_before_);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(b_spline_after);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(BSpline1DFig5_26, RemovesKnot1_0CorrectlyOneTime) {  // NOLINT
  ASSERT_THAT(bspline_1d_after_->RemoveKnot(ParametricCoordinate(1), 0, 0.1), 1);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 1);
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[6].Get(), DoubleEq(2));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() - 1);
  std::vector<spl::ControlPoint> new_control_points = {
      spl::ControlPoint({0.0, 0.0}), spl::ControlPoint({0.0, 1.5}), spl::ControlPoint({1.0, 2.0}),
      spl::ControlPoint({3.0, 2.0}), spl::ControlPoint({4.0, 1.5}), spl::ControlPoint({4.0, 0.0})
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    ASSERT_THAT(bspline_1d_after_->GetControlPoint(i, 0),
                DoubleEq(new_control_points[i].GetValueForDimension(Dimension{0})));
    ASSERT_THAT(bspline_1d_after_->GetControlPoint(i, 1),
                DoubleEq(new_control_points[i].GetValueForDimension(Dimension{1})));
  }
  ASSERT_THAT(bspline_1d_before_->AreGeometricallyEqual(*bspline_1d_after_), true);
}

TEST_F(BSpline1DFig5_26, RemovesKnot1_0CorrectlyTwoTimes) {  // NOLINT
  ASSERT_THAT(bspline_1d_after_->RemoveKnot(ParametricCoordinate(1), 0, 0.0, 2), 2);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 2);
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[5].Get(), DoubleEq(2));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() - 2);
  std::vector<spl::ControlPoint> new_control_points = {
      spl::ControlPoint({0.0, 0.0}), spl::ControlPoint({0.0, 1.5}), spl::ControlPoint({2.0, 2.5}),
      spl::ControlPoint({4.0, 1.5}), spl::ControlPoint({4.0, 0.0})
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    ASSERT_THAT(bspline_1d_after_->GetControlPoint(i, 0),
                DoubleEq(new_control_points[i].GetValueForDimension(Dimension{0})));
    ASSERT_THAT(bspline_1d_after_->GetControlPoint(i, 1),
                DoubleEq(new_control_points[i].GetValueForDimension(Dimension{1})));
  }
  ASSERT_THAT(bspline_1d_before_->AreGeometricallyEqual(*bspline_1d_after_), true);
}

TEST_F(BSpline1DFig5_26, RemovesKnot1_0CorrectlyThreeTimesAtOnce) {  // NOLINT
  ASSERT_THAT(bspline_1d_after_->RemoveKnot(ParametricCoordinate(1), 0, 0.5, 3), 3);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 3);
  ASSERT_THAT((*bspline_1d_after_->GetKnotVector(0))[4].Get(), DoubleEq(2));
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() - 3);
  std::vector<spl::ControlPoint> new_control_points = {
      spl::ControlPoint({0.0, 0.0}), spl::ControlPoint({0.0, 3.0}), spl::ControlPoint({4.0, 3.0}),
      spl::ControlPoint({4.0, 0.0})
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    ASSERT_THAT(bspline_1d_after_->GetControlPoint(i, 0),
                DoubleEq(new_control_points[i].GetValueForDimension(Dimension{0})));
    ASSERT_THAT(bspline_1d_after_->GetControlPoint(i, 1),
                DoubleEq(new_control_points[i].GetValueForDimension(Dimension{1})));
  }
  ASSERT_THAT(bspline_1d_before_->AreGeometricallyEqual(*bspline_1d_after_, 0.2), false);
  ASSERT_THAT(bspline_1d_before_->AreGeometricallyEqual(*bspline_1d_after_, 0.3), true);
}

TEST_F(BSpline1DFig5_26, RemovesOnlyTwoKnots1_0CorrectlyWithTolerance0_1) {  // NOLINT
  ASSERT_THAT(bspline_1d_after_->RemoveKnot(ParametricCoordinate(1), 0, 0.1, 3), 2);
}

class NURBS1DFig5_26 : public Test {  // NOLINT
 public:
  NURBS1DFig5_26() {
    std::array<Degree, 1> degree = {Degree{3}};
    baf::KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(baf::KnotVector(
        {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0},
         ParametricCoordinate{1},
         ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{2}, ParametricCoordinate{2},
         ParametricCoordinate{2},
         ParametricCoordinate{2}}))};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint(std::vector<double>({0.0, 0.0})),
        spl::ControlPoint(std::vector<double>({0.0, 1.5})),
        spl::ControlPoint(std::vector<double>({1.0, 2.0})),
        spl::ControlPoint(std::vector<double>({2.0, 2.0})),
        spl::ControlPoint(std::vector<double>({3.0, 2.0})),
        spl::ControlPoint(std::vector<double>({4.0, 1.5})),
        spl::ControlPoint(std::vector<double>({4.0, 0.0}))
    };
    std::vector<double> weights = {1, 0.8, 1.4, 1.2, 1.1, 0.8, 1.4};
    nurbs_1d_before_ = std::make_shared<spl::NURBS<1>>(knot_vector_before, degree, control_points, weights);
    spl::NURBS<1> nurbs_after(*nurbs_1d_before_);
    nurbs_1d_after_ = std::make_shared<spl::NURBS<1>>(nurbs_after);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_before_;
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_after_;
};

TEST_F(NURBS1DFig5_26, RemovesKnot1_0CorrectlyOneTime) {  // NOLINT
  nurbs_1d_after_->RemoveKnot(ParametricCoordinate(1), 0, 0.81);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_1d_before_->GetKnotVector(0)->GetNumberOfKnots() - 1);
  ASSERT_THAT((*nurbs_1d_after_->GetKnotVector(0))[6].Get(), DoubleEq(2));
  ASSERT_THAT(nurbs_1d_after_->GetNumberOfControlPoints(), nurbs_1d_before_->GetNumberOfControlPoints() - 1);
  std::vector<spl::ControlPoint> new_control_points = {
      spl::ControlPoint(std::vector<double>({0.0, 0.0})),
      spl::ControlPoint(std::vector<double>({0.0, 1.5})),
      spl::ControlPoint(std::vector<double>({1.0, 2.0})),
      spl::ControlPoint(std::vector<double>({3.0, 2.0})),
      spl::ControlPoint(std::vector<double>({4.0, 1.5})),
      spl::ControlPoint(std::vector<double>({4.0, 0.0}))
  };
  std::vector<double> new_weights = {1, 0.8, 1.4, 1.1, 0.8, 1.4};
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(nurbs_1d_after_->GetControlPoint(i, j),
                  DoubleEq(new_control_points[i].GetValueForDimension(Dimension{j})));
    }
    ASSERT_THAT(nurbs_1d_after_->GetWeight(i), DoubleEq(new_weights[i]));
  }
  ASSERT_THAT(nurbs_1d_before_->AreGeometricallyEqual(*nurbs_1d_after_, 0.1), false);
  ASSERT_THAT(nurbs_1d_before_->AreGeometricallyEqual(*nurbs_1d_after_, 0.13), true);
}

TEST_F(NURBS1DFig5_26, RemovesOnlyOneKnot1_0CorrectlyWithTolerance0_81) {  // NOLINT
  ASSERT_THAT(nurbs_1d_after_->RemoveKnot(ParametricCoordinate(1), 0, 0.81, 2), 1);
}

class BSpline2DFig5_28 : public Test {  // NOLINT
 public:
  BSpline2DFig5_28() {
    std::array<Degree, 2> degree = {Degree{2}, Degree{3}};
    ParametricCoordinate zero(0), one(1);
    baf::KnotVectors<2> knot_vector_before = {
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {zero, zero, zero, ParametricCoordinate{0.25}, ParametricCoordinate{0.5}, ParametricCoordinate{0.75}, one,
             one, one})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {zero, zero, zero, zero, ParametricCoordinate{0.3}, ParametricCoordinate{0.3}, ParametricCoordinate{0.3},
             ParametricCoordinate{0.7},
             one, one, one, one}))};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint({0, 1, 4}), spl::ControlPoint({1, 1, 4.5}), spl::ControlPoint({2.5, 1, 3.5}),
        spl::ControlPoint({3.5, 1, 2.5}), spl::ControlPoint({5, 1, 2}), spl::ControlPoint({6.5, 1, 2.5}),

        spl::ControlPoint({0, 1.5, 3.9}), spl::ControlPoint({1, 1.5, 4.4}), spl::ControlPoint({2.5, 1.5, 3.4}),
        spl::ControlPoint({3.5, 1.5, 2.4}), spl::ControlPoint({5, 1.5, 1.9}), spl::ControlPoint({6.5, 1.5, 2.4}),

        spl::ControlPoint({0, 2, 3.8}), spl::ControlPoint({1, 2, 4.3}), spl::ControlPoint({2.5, 2, 3.3}),
        spl::ControlPoint({3.5, 2, 2.3}), spl::ControlPoint({5, 2, 1.8}), spl::ControlPoint({6.5, 2, 2.3}),

        spl::ControlPoint({0, 2.5, 3.7}), spl::ControlPoint({1, 2.5, 4.2}), spl::ControlPoint({2.5, 2.5, 3.2}),
        spl::ControlPoint({3.5, 2.5, 2.2}), spl::ControlPoint({5, 2.5, 1.7}), spl::ControlPoint({6.5, 2.5, 2.2}),

        spl::ControlPoint({0, 3, 3.6}), spl::ControlPoint({1, 3, 4.1}), spl::ControlPoint({2.5, 3, 3.1}),
        spl::ControlPoint({3.5, 3, 2.1}), spl::ControlPoint({5, 3, 1.6}), spl::ControlPoint({6.5, 3, 2.1}),

        spl::ControlPoint({0, 4, 3.3}), spl::ControlPoint({1, 4, 3.8}), spl::ControlPoint({2.5, 4, 2.8}),
        spl::ControlPoint({3.5, 4, 1.8}), spl::ControlPoint({5, 4, 1.3}), spl::ControlPoint({6.5, 4, 1.8}),

        spl::ControlPoint({0, 5, 3}), spl::ControlPoint({1, 5, 3.5}), spl::ControlPoint({2.5, 5, 2.5}),
        spl::ControlPoint({3.5, 5, 1.5}), spl::ControlPoint({5, 5, 1}), spl::ControlPoint({6.5, 5, 1.5}),

        spl::ControlPoint({0, 5.5, 3}), spl::ControlPoint({1, 5.5, 3.5}), spl::ControlPoint({2.5, 5.5, 2.5}),
        spl::ControlPoint({3.5, 5.5, 1.5}), spl::ControlPoint({5, 5.5, 1}), spl::ControlPoint({6.5, 5.5, 1.5})
    };
    bspline_2d_before_ = std::make_shared<spl::BSpline<2>>(knot_vector_before, degree, control_points);
    spl::BSpline<2> b_spline_after(*bspline_2d_before_);
    bspline_2d_after_ = std::make_shared<spl::BSpline<2>>(b_spline_after);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> bspline_2d_before_;
  std::shared_ptr<spl::BSpline<2>> bspline_2d_after_;
};

TEST_F(BSpline2DFig5_28, RemovesKnot0_3CorrectlyOneTime) {  // NOLINT
  ASSERT_THAT(bspline_2d_after_->RemoveKnot(ParametricCoordinate(0.3), 1, 0.075), 1);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(1)->GetNumberOfKnots() - 1);
  ASSERT_THAT((*bspline_2d_after_->GetKnotVector(1))[6].Get(), DoubleEq(0.7));
  ASSERT_THAT(bspline_2d_after_->GetNumberOfControlPoints(), bspline_2d_before_->GetNumberOfControlPoints() - 6);
  ASSERT_THAT(bspline_2d_before_->AreGeometricallyEqual(*bspline_2d_after_), false);
  ASSERT_THAT(bspline_2d_before_->AreGeometricallyEqual(*bspline_2d_after_, 0.075), true);
}

TEST_F(BSpline2DFig5_28, RemovesKnot0_3CorrectlyTwoTimes) {  // NOLINT
  ASSERT_THAT(bspline_2d_after_->RemoveKnot(ParametricCoordinate(0.3), 1, 0.12, 2), 2);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(1)->GetNumberOfKnots() - 2);
  ASSERT_THAT((*bspline_2d_after_->GetKnotVector(1))[5].Get(), DoubleEq(0.7));
  ASSERT_THAT(bspline_2d_after_->GetNumberOfControlPoints(), bspline_2d_before_->GetNumberOfControlPoints() - 12);
  ASSERT_THAT(bspline_2d_before_->AreGeometricallyEqual(*bspline_2d_after_), false);
  ASSERT_THAT(bspline_2d_before_->AreGeometricallyEqual(*bspline_2d_after_, 0.1), true);
}

class NURBS2DFig5_28 : public Test {  // NOLINT
 public:
  NURBS2DFig5_28() {
    std::array<Degree, 2> degree = {Degree{2}, Degree{3}};
    ParametricCoordinate zero(0), one(1);
    baf::KnotVectors<2> knot_vector_before = {
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {zero, zero, zero, ParametricCoordinate{0.25}, ParametricCoordinate{0.5}, ParametricCoordinate{0.75}, one,
             one, one})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {zero, zero, zero, zero, ParametricCoordinate{0.3}, ParametricCoordinate{0.3}, ParametricCoordinate{0.3},
             ParametricCoordinate{0.7},
             one, one, one, one}))};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint({0, 1, 4}), spl::ControlPoint({1, 1, 4.5}), spl::ControlPoint({2.5, 1, 3.5}),
        spl::ControlPoint({3.5, 1, 2.5}), spl::ControlPoint({5, 1, 2}), spl::ControlPoint({6.5, 1, 2.5}),

        spl::ControlPoint({0, 1.5, 3.9}), spl::ControlPoint({1, 1.5, 4.4}), spl::ControlPoint({2.5, 1.5, 3.4}),
        spl::ControlPoint({3.5, 1.5, 2.4}), spl::ControlPoint({5, 1.5, 1.9}), spl::ControlPoint({6.5, 1.5, 2.4}),

        spl::ControlPoint({0, 2, 3.8}), spl::ControlPoint({1, 2, 4.3}), spl::ControlPoint({2.5, 2, 3.3}),
        spl::ControlPoint({3.5, 2, 2.3}), spl::ControlPoint({5, 2, 1.8}), spl::ControlPoint({6.5, 2, 2.3}),

        spl::ControlPoint({0, 2.5, 3.7}), spl::ControlPoint({1, 2.5, 4.2}), spl::ControlPoint({2.5, 2.5, 3.2}),
        spl::ControlPoint({3.5, 2.5, 2.2}), spl::ControlPoint({5, 2.5, 1.7}), spl::ControlPoint({6.5, 2.5, 2.2}),

        spl::ControlPoint({0, 3, 3.6}), spl::ControlPoint({1, 3, 4.1}), spl::ControlPoint({2.5, 3, 3.1}),
        spl::ControlPoint({3.5, 3, 2.1}), spl::ControlPoint({5, 3, 1.6}), spl::ControlPoint({6.5, 3, 2.1}),

        spl::ControlPoint({0, 4, 3.3}), spl::ControlPoint({1, 4, 3.8}), spl::ControlPoint({2.5, 4, 2.8}),
        spl::ControlPoint({3.5, 4, 1.8}), spl::ControlPoint({5, 4, 1.3}), spl::ControlPoint({6.5, 4, 1.8}),

        spl::ControlPoint({0, 5, 3}), spl::ControlPoint({1, 5, 3.5}), spl::ControlPoint({2.5, 5, 2.5}),
        spl::ControlPoint({3.5, 5, 1.5}), spl::ControlPoint({5, 5, 1}), spl::ControlPoint({6.5, 5, 1.5}),

        spl::ControlPoint({0, 5.5, 3}), spl::ControlPoint({1, 5.5, 3.5}), spl::ControlPoint({2.5, 5.5, 2.5}),
        spl::ControlPoint({3.5, 5.5, 1.5}), spl::ControlPoint({5, 5.5, 1}), spl::ControlPoint({6.5, 5.5, 1.5})
    };
    std::vector<double> weights = {1, 0.98, 0.96, 0.94, 0.96, 0.94,
                                   1, 0.99, 0.97, 0.96, 0.98, 1,
                                   1.01, 1.02, 1.01, 1, 1, 1,
                                   1, 1, 1, 1, 1.01, 1.02,
                                   1.01, 1.03, 1.04, 1.02, 1.01, 1,
                                   1, 0.95, 0.98, 1, 1.05, 1.01,
                                   0.95, 0.9, 0.8, 0.9, 1, 0.95,
                                   0.9, 1, 1.01, 1.02, 1, 0.8};
    nurbs_2d_before_ = std::make_shared<spl::NURBS<2>>(knot_vector_before, degree, control_points, weights);
    spl::NURBS<2> nurbs_after(*nurbs_2d_before_);
    nurbs_2d_after_ = std::make_shared<spl::NURBS<2>>(nurbs_after);
  }

 protected:
  std::shared_ptr<spl::NURBS<2>> nurbs_2d_before_;
  std::shared_ptr<spl::NURBS<2>> nurbs_2d_after_;
};

TEST_F(NURBS2DFig5_28, RemovesKnot0_3CorrectlyOneTime) {  // NOLINT
  ASSERT_THAT(nurbs_2d_after_->RemoveKnot(ParametricCoordinate(0.3), 1, 2.2), 1);
  ASSERT_THAT(nurbs_2d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              nurbs_2d_before_->GetKnotVector(1)->GetNumberOfKnots() - 1);
  ASSERT_THAT((*nurbs_2d_after_->GetKnotVector(1))[6].Get(), DoubleEq(0.7));
  ASSERT_THAT(nurbs_2d_after_->GetNumberOfControlPoints(), nurbs_2d_before_->GetNumberOfControlPoints() - 6);
  ASSERT_THAT(nurbs_2d_before_->AreGeometricallyEqual(*nurbs_2d_after_), false);
  ASSERT_THAT(nurbs_2d_before_->AreGeometricallyEqual(*nurbs_2d_after_, 0.08), true);
}

TEST_F(NURBS2DFig5_28, RemovesKnot0_3CorrectlyTwoTimes) {  // NOLINT
  ASSERT_THAT(nurbs_2d_after_->RemoveKnot(ParametricCoordinate(0.3), 1, 2.2, 2), 2);
  ASSERT_THAT(nurbs_2d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              nurbs_2d_before_->GetKnotVector(1)->GetNumberOfKnots() - 2);
  ASSERT_THAT((*nurbs_2d_after_->GetKnotVector(1))[5].Get(), DoubleEq(0.7));
  ASSERT_THAT(nurbs_2d_after_->GetNumberOfControlPoints(), nurbs_2d_before_->GetNumberOfControlPoints() - 12);
  ASSERT_THAT(nurbs_2d_before_->AreGeometricallyEqual(*nurbs_2d_after_), false);
  ASSERT_THAT(nurbs_2d_before_->AreGeometricallyEqual(*nurbs_2d_after_, 0.07), true);
}

class A3DBSplineForKnotRemoval : public Test {  // NOLINT
 public:
  A3DBSplineForKnotRemoval() {
    std::array<Degree, 3> degree = {Degree{2}, Degree{1}, Degree{2}};
    ParametricCoordinate zero(0), one(1);
    baf::KnotVectors<3> knot_vector_before = {
        std::make_shared<baf::KnotVector>(baf::KnotVector({zero, zero, zero, one, one, one})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({zero, zero, one, one})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({zero, zero, zero, ParametricCoordinate{0.3}, ParametricCoordinate{0.3}, one, one, one}))};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint({0, 1, 4}), spl::ControlPoint({1, 1, 4.5}), spl::ControlPoint({2.5, 1, 3.5}),
        spl::ControlPoint({0, 1, 6}), spl::ControlPoint({1, 1, 6.5}), spl::ControlPoint({2.5, 1, 5.5}),

        spl::ControlPoint({0, 1.5, 3.9}), spl::ControlPoint({1, 1.5, 4.4}), spl::ControlPoint({2.5, 1.5, 3.4}),
        spl::ControlPoint({0, 1.5, 5.9}), spl::ControlPoint({1, 1.5, 6.4}), spl::ControlPoint({2.5, 1.5, 5.4}),

        spl::ControlPoint({0, 2, 3.8}), spl::ControlPoint({1, 2, 4.3}), spl::ControlPoint({2.5, 2, 3.3}),
        spl::ControlPoint({0, 2, 5.8}), spl::ControlPoint({1, 2, 6.3}), spl::ControlPoint({2.5, 2, 5.3}),

        spl::ControlPoint({0, 2.5, 3.7}), spl::ControlPoint({1, 2.5, 4.2}), spl::ControlPoint({2.5, 2.5, 3.2}),
        spl::ControlPoint({0, 2.5, 5.7}), spl::ControlPoint({1, 2.5, 6.2}), spl::ControlPoint({2.5, 2.5, 5.2}),

        spl::ControlPoint({0, 3, 3.6}), spl::ControlPoint({1, 3, 4.1}), spl::ControlPoint({2.5, 3, 3.1}),
        spl::ControlPoint({0, 3, 5.6}), spl::ControlPoint({1, 3, 6.1}), spl::ControlPoint({2.5, 3, 5.1})
    };
    bspline_3d_before_ = std::make_shared<spl::BSpline<3>>(knot_vector_before, degree, control_points);
    spl::BSpline<3> b_spline_after(*bspline_3d_before_);
    bspline_3d_after_ = std::make_shared<spl::BSpline<3>>(b_spline_after);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> bspline_3d_before_;
  std::shared_ptr<spl::BSpline<3>> bspline_3d_after_;
};

TEST_F(A3DBSplineForKnotRemoval, RemovesKnot0_3CorrectlyOneTime) {  // NOLINT
  ASSERT_THAT(bspline_3d_after_->RemoveKnot(ParametricCoordinate(0.3), 2, 0.21), 1);
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(2)->GetNumberOfKnots(),
              bspline_3d_before_->GetKnotVector(2)->GetNumberOfKnots() - 1);
  ASSERT_THAT((*bspline_3d_after_->GetKnotVector(2))[5].Get(), DoubleEq(1));
  ASSERT_THAT(bspline_3d_after_->GetNumberOfControlPoints(), bspline_3d_before_->GetNumberOfControlPoints() - 6);
  ASSERT_THAT(bspline_3d_before_->AreGeometricallyEqual(*bspline_3d_after_, 0.1), false);
  ASSERT_THAT(bspline_3d_before_->AreGeometricallyEqual(*bspline_3d_after_, 0.2), true);
}

TEST_F(A3DBSplineForKnotRemoval, RemovesKnot0_3CorrectlyTwoTimes) {  // NOLINT
  ASSERT_THAT(bspline_3d_after_->RemoveKnot(ParametricCoordinate(0.3), 2, 0.28, 2), 2);
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(2)->GetNumberOfKnots(),
              bspline_3d_before_->GetKnotVector(2)->GetNumberOfKnots() - 2);
  ASSERT_THAT((*bspline_3d_after_->GetKnotVector(2))[4].Get(), DoubleEq(1));
  ASSERT_THAT(bspline_3d_after_->GetNumberOfControlPoints(), bspline_3d_before_->GetNumberOfControlPoints() - 12);
  ASSERT_THAT(bspline_3d_before_->AreGeometricallyEqual(*bspline_3d_after_, 0.1), false);
  ASSERT_THAT(bspline_3d_before_->AreGeometricallyEqual(*bspline_3d_after_, 0.2), true);
}

class A3DNURBSForKnotRemoval : public Test {  // NOLINT
 public:
  A3DNURBSForKnotRemoval() {
    std::array<Degree, 3> degree = {Degree{2}, Degree{1}, Degree{2}};
    ParametricCoordinate zero(0), one(1);
    baf::KnotVectors<3> knot_vector_before = {
        std::make_shared<baf::KnotVector>(baf::KnotVector({zero, zero, zero, one, one, one})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({zero, zero, one, one})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({zero, zero, zero, ParametricCoordinate{0.3}, ParametricCoordinate{0.3}, one, one, one}))};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint({0, 1, 4}), spl::ControlPoint({1, 1, 4.5}), spl::ControlPoint({2.5, 1, 3.5}),
        spl::ControlPoint({0, 1, 6}), spl::ControlPoint({1, 1, 6.5}), spl::ControlPoint({2.5, 1, 5.5}),

        spl::ControlPoint({0, 1.5, 3.9}), spl::ControlPoint({1, 1.5, 4.4}), spl::ControlPoint({2.5, 1.5, 3.4}),
        spl::ControlPoint({0, 1.5, 5.9}), spl::ControlPoint({1, 1.5, 6.4}), spl::ControlPoint({2.5, 1.5, 5.4}),

        spl::ControlPoint({0, 2, 3.8}), spl::ControlPoint({1, 2, 4.3}), spl::ControlPoint({2.5, 2, 3.3}),
        spl::ControlPoint({0, 2, 5.8}), spl::ControlPoint({1, 2, 6.3}), spl::ControlPoint({2.5, 2, 5.3}),

        spl::ControlPoint({0, 2.5, 3.7}), spl::ControlPoint({1, 2.5, 4.2}), spl::ControlPoint({2.5, 2.5, 3.2}),
        spl::ControlPoint({0, 2.5, 5.7}), spl::ControlPoint({1, 2.5, 6.2}), spl::ControlPoint({2.5, 2.5, 5.2}),

        spl::ControlPoint({0, 3, 3.6}), spl::ControlPoint({1, 3, 4.1}), spl::ControlPoint({2.5, 3, 3.1}),
        spl::ControlPoint({0, 3, 5.6}), spl::ControlPoint({1, 3, 6.1}), spl::ControlPoint({2.5, 3, 5.1})
    };
    std::vector<double> weights = {1, 1.02, 1.04, 1.06, 1.04, 1.02,
                                   0.98, 0.99, 1, 0.99, 1, 1.01,
                                   1, 1.02, 0.99, 1.05, 0.98, 1,
                                   0.95, 0.96, 0.97, 0.98, 0.99, 1,
                                   0.98, 1, 1.03, 1.05, 1, 0.97};
    nurbs_3d_before_ = std::make_shared<spl::NURBS<3>>(knot_vector_before, degree, control_points, weights);
    spl::NURBS<3> nurbs_after(*nurbs_3d_before_);
    nurbs_3d_after_ = std::make_shared<spl::NURBS<3>>(nurbs_after);
  }

 protected:
  std::shared_ptr<spl::NURBS<3>> nurbs_3d_before_;
  std::shared_ptr<spl::NURBS<3>> nurbs_3d_after_;
};

TEST_F(A3DNURBSForKnotRemoval, RemovesKnot0_3CorrectlyOneTime) {  // NOLINT
  ASSERT_THAT(nurbs_3d_after_->RemoveKnot(ParametricCoordinate(0.3), 2, 3.86), 1);
  ASSERT_THAT(nurbs_3d_after_->GetKnotVector(2)->GetNumberOfKnots(),
              nurbs_3d_before_->GetKnotVector(2)->GetNumberOfKnots() - 1);
  ASSERT_THAT((*nurbs_3d_after_->GetKnotVector(2))[5].Get(), DoubleEq(1));
  ASSERT_THAT(nurbs_3d_after_->GetNumberOfControlPoints(), nurbs_3d_before_->GetNumberOfControlPoints() - 6);
  ASSERT_THAT(nurbs_3d_before_->AreGeometricallyEqual(*nurbs_3d_after_, 0.1), false);
  ASSERT_THAT(nurbs_3d_before_->AreGeometricallyEqual(*nurbs_3d_after_, 0.2), true);
}

TEST_F(A3DNURBSForKnotRemoval, RemovesKnot0_3CorrectlyTwoTimes) {  // NOLINT
  ASSERT_THAT(nurbs_3d_after_->RemoveKnot(ParametricCoordinate(0.3), 2, 5.3, 2), 2);
  ASSERT_THAT(nurbs_3d_after_->GetKnotVector(2)->GetNumberOfKnots(),
              nurbs_3d_before_->GetKnotVector(2)->GetNumberOfKnots() - 2);
  ASSERT_THAT((*nurbs_3d_after_->GetKnotVector(2))[4].Get(), DoubleEq(1));
  ASSERT_THAT(nurbs_3d_after_->GetNumberOfControlPoints(), nurbs_3d_before_->GetNumberOfControlPoints() - 12);
  ASSERT_THAT(nurbs_3d_before_->AreGeometricallyEqual(*nurbs_3d_after_, 0.1), false);
  ASSERT_THAT(nurbs_3d_before_->AreGeometricallyEqual(*nurbs_3d_after_, 0.22), true);
}
