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

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

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
