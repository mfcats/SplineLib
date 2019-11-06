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

#include "src/spl/b_spline.h"
#include "src/spl/nurbs.h"
#include "src/util/any_casts.h"

using testing::Test;
using testing::DoubleEq;

using namespace splinelib::src;

class AnyCasts : public Test {
 public:
  AnyCasts() = default;
};

TEST_F(AnyCasts, ThrowLogicErrorForIntegerAsInputInsteadOfSplinePointer) {  // NOLINT
  ASSERT_THROW(util::any_casts::GetSplineDimension<4>(std::make_any<int>(8)), std::logic_error);
}

TEST_F(AnyCasts, ThrowInvalidArgumentForGetSplineWithIntegerAsArgument) {  // NOLINT
  ASSERT_THROW(util::any_casts::GetSpline<1>(std::make_any<int>(8)), std::invalid_argument);
}

TEST_F(AnyCasts, ThrowInvalidArgumentForIsRationalWithIntegerAsArgument) {  // NOLINT
  ASSERT_THROW(util::any_casts::IsRational<1>(std::make_any<int>(8)), std::invalid_argument);
}

class A1DAnyBSplineForAnyCasts : public Test {
 public:
  A1DAnyBSplineForAnyCasts() {
    std::array<Degree, 1> degree = {Degree(1)};
    baf::KnotVectors<1> knot_vector_ptr = {std::make_shared<baf::KnotVector>(baf::KnotVector(
        {ParametricCoordinate(0), ParametricCoordinate(0), ParametricCoordinate(1), ParametricCoordinate(1)}))};
    std::vector<baf::ControlPoint> control_points = {baf::ControlPoint{0.0}, baf::ControlPoint{1.2}};
    spl::BSpline<1> b_spline_1d(knot_vector_ptr, degree, control_points);
    std::shared_ptr<spl::BSpline<1>> b_spline_1d_ptr = std::make_shared<spl::BSpline<1>>(b_spline_1d);
    b_spline_1d_any_ = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_1d_ptr);
  }

 protected:
  std::any b_spline_1d_any_;
};

TEST_F(A1DAnyBSplineForAnyCasts, Returns1_2AsSecondControlPointAfterCastTo1DBSpline) {  // NOLINT
  ASSERT_THAT(util::any_casts::GetSpline<1>(b_spline_1d_any_)->GetControlPoint(1).GetValue(0), DoubleEq(1.2));
}

TEST_F(A1DAnyBSplineForAnyCasts, ReturnsSplineDimension1) {  // NOLINT
  ASSERT_THAT(util::any_casts::GetSplineDimension<1>(b_spline_1d_any_), 1);
}

TEST_F(A1DAnyBSplineForAnyCasts, ReturnsThatTheSplineIsNotRational) {  // NOLINT
  ASSERT_THAT(util::any_casts::IsRational<1>(b_spline_1d_any_), false);
}

class A2DAnyNURBSForAnyCasts : public Test {
 public:
  A2DAnyNURBSForAnyCasts() {
    std::array<Degree, 2> degree = {Degree(1), Degree(1)};
    baf::KnotVectors<2> knot_vector_ptr = {
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParametricCoordinate(0), ParametricCoordinate(0),
                                                           ParametricCoordinate(1), ParametricCoordinate(1)})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParametricCoordinate(0), ParametricCoordinate(0),
                                                           ParametricCoordinate(1), ParametricCoordinate(1)}))};
    std::vector<baf::ControlPoint> control_points = {baf::ControlPoint{0.0}, baf::ControlPoint{1.2},
                                                     baf::ControlPoint{-1.0}, baf::ControlPoint{-2.2}};
    std::vector<double> weights = {1.0, 2.0, 1.5, 0.7};
    spl::NURBS<2> nurbs_2d(knot_vector_ptr, degree, control_points, weights);
    std::shared_ptr<spl::NURBS<2>> nurbs_2d_ptr = std::make_shared<spl::NURBS<2>>(nurbs_2d);
    nurbs_2d_any_ = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbs_2d_ptr);
  }

 protected:
  std::any nurbs_2d_any_;
};

TEST_F(A2DAnyNURBSForAnyCasts, ReturnsMinus1AsThirdControlPointAfterCastTo2DNURBS) {  // NOLINT
  ASSERT_THAT(util::any_casts::GetSpline<2>(nurbs_2d_any_)->GetControlPoint(2).GetValue(0), DoubleEq(-1.0));
}

TEST_F(A2DAnyNURBSForAnyCasts, ReturnsSplineDimension2) {  // NOLINT
  ASSERT_THAT(util::any_casts::GetSplineDimension<2>(nurbs_2d_any_), 2);
}

TEST_F(A2DAnyNURBSForAnyCasts, ReturnsThatTheSplineIsRational) {  // NOLINT
  ASSERT_THAT(util::any_casts::IsRational<2>(nurbs_2d_any_), true);
}
