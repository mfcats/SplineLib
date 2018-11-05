/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <nurbs.h>
#include "gmock/gmock.h"

#include "b_spline.h"

using testing::Test;
using testing::DoubleEq;

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
    nurbs_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    nurbs_1d_after_ = std::make_shared<spl::BSpline<1>>(knot_vector_after, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> nurbs_1d_before_;
  std::shared_ptr<spl::BSpline<1>> nurbs_1d_after_;
};

TEST_F(BSplineEx5_1, InsertsKnot2_5Correctly) {  // NOLINT
  nurbs_1d_after_->InsertKnot(ParamCoord(2.5), 0);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetKnot(6).get(), DoubleEq(2.5));
  ASSERT_THAT(nurbs_1d_after_->GetNumberOfControlPoints(), nurbs_1d_before_->GetNumberOfControlPoints() + 1);
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
      ASSERT_THAT(nurbs_1d_after_->GetControlPoint({i}, j), DoubleEq(new_control_points[i].GetValue(j)));
    }
  }
  for (int i = 0; i <= 50; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / 10.0)};
    ASSERT_THAT(nurbs_1d_after_->Evaluate(param_coord, {0})[0],
                DoubleEq(nurbs_1d_before_->Evaluate(param_coord, {0})[0]));
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
  std::vector<baf::ControlPoint> new_control_points = {
      baf::ControlPoint(std::vector<double>({0.0, 1.0})),
      baf::ControlPoint(std::vector<double>({1.0, 2.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({7.0 / 3.0, 2.0 / 3.0})),
      baf::ControlPoint(std::vector<double>({5.0, 2.0})),
      baf::ControlPoint(std::vector<double>({1.0, 1.0})),
      baf::ControlPoint(std::vector<double>({1.0, 2.0})),
      baf::ControlPoint(std::vector<double>({1.3, 2.3}))
  };
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(nurbs_1d_after_->GetControlPoint({i}, j), DoubleEq(new_control_points[i].GetValue(j)));
    }
  }
  ASSERT_THAT(nurbs_1d_after_->GetWeights().size(), nurbs_1d_before_->GetWeights().size() + 1);
  std::vector<double> new_weights = {1.0, 1.0, 2.0, 10.0 / 3.0, 3.0, 1.0, 4.0, 2.0, 1.0};
  for (int i = 0; i < static_cast<int>(new_control_points.size()); ++i) {
    ASSERT_THAT(nurbs_1d_after_->GetWeight({i}), DoubleEq(new_weights[i]));
  }
//  for (int i = 0; i <= 50; ++i) {
//    std::array<ParamCoord, 1> param_coord{ParamCoord(i / 10.0)};
//    ASSERT_THAT(nurbs_1d_after_->Evaluate(param_coord, {0})[0],
//                DoubleEq(nurbs_1d_before_->Evaluate(param_coord, {0})[0]));
//  }
}

