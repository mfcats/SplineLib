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

using testing::Test;
using testing::DoubleEq;

class A1DBSplineForKnotInsertion : public Test {  // NOLINT
 public:
  A1DBSplineForKnotInsertion() {
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

TEST_F(A1DBSplineForKnotInsertion, InsertsKnot2_5Correctly) {  // NOLINT
  ASSERT_THAT(bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots(), 12);
  bspline_1d_after_->InsertKnot(ParamCoord(2.5), 0);
  ASSERT_THAT(bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots(), 12);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() + 1);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetKnot(6).get(), DoubleEq(2.5));
  for (int i = 0; i <= 50; ++i) {
    ASSERT_THAT(bspline_1d_after_->Evaluate({ParamCoord(i / 10.0)}, {0})[0],
                DoubleEq(bspline_1d_before_->Evaluate({ParamCoord(i / 10.0)}, {0})[0]));
  }
}
