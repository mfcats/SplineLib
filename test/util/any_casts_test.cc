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
#include "src/util/any_casts.h"

using testing::DoubleEq;
using testing::Test;

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
    std::shared_ptr<spl::BSpline<1>> b_spline_1d_ptr = std::make_shared<spl::BSpline<1>>(knot_vector_ptr, degree,
                                                                                         control_points);
    b_spline_1d_any_ = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_1d_ptr);
  }

 protected:
  std::any b_spline_1d_any_;
};

TEST_F(A1DAnyBSplineForAnyCasts, Returns1_2AsSecondControlPointAfterCastTo1DBSpline) {  // NOLINT
  ASSERT_THAT(util::any_casts::GetSpline<1>(b_spline_1d_any_)->GetControlPoint(1).GetValueForDimension(Dimension{0}),
      DoubleEq(1.2));
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
    baf::KnotVector knot_vector({ParametricCoordinate(0), ParametricCoordinate(0), ParametricCoordinate(1),
                                 ParametricCoordinate(1)});
    baf::KnotVectors<2> knot_vector_ptr = {std::make_shared<baf::KnotVector>(knot_vector),
                                           std::make_shared<baf::KnotVector>(knot_vector)};
    std::vector<baf::ControlPoint> control_points = {baf::ControlPoint{0.0}, baf::ControlPoint{1.2},
                                                     baf::ControlPoint{-1.0}, baf::ControlPoint{-2.2}};
    std::vector<double> weights = {1.0, 2.0, 1.5, 0.7};
    std::shared_ptr<spl::NURBS<2>> nurbs_2d_ptr = std::make_shared<spl::NURBS<2>>(knot_vector_ptr, degree,
                                                                                  control_points, weights);
    nurbs_2d_any_ = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbs_2d_ptr);
  }

 protected:
  std::any nurbs_2d_any_;
};

TEST_F(A2DAnyNURBSForAnyCasts, ReturnsMinus1AsThirdControlPointAfterCastTo2DNURBS) {  // NOLINT
  ASSERT_THAT(util::any_casts::GetSpline<2>(nurbs_2d_any_)->GetControlPoint(2).GetValueForDimension(Dimension{0}),
      DoubleEq(-1.0));
}

TEST_F(A2DAnyNURBSForAnyCasts, ReturnsSplineDimension2) {  // NOLINT
  ASSERT_THAT(util::any_casts::GetSplineDimension<2>(nurbs_2d_any_), 2);
}

TEST_F(A2DAnyNURBSForAnyCasts, ReturnsThatTheSplineIsRational) {  // NOLINT
  ASSERT_THAT(util::any_casts::IsRational<2>(nurbs_2d_any_), true);
}

class A3DAnyBSplineForAnyCasts : public Test {
 public:
  A3DAnyBSplineForAnyCasts() {
    std::array<Degree, 3> degree = {Degree(1), Degree(1), Degree(1)};
    baf::KnotVector knot_vector({ParametricCoordinate(0), ParametricCoordinate(0), ParametricCoordinate(1),
                                 ParametricCoordinate(1)});
    baf::KnotVectors<3> knot_vector_ptr = {std::make_shared<baf::KnotVector>(knot_vector),
                                           std::make_shared<baf::KnotVector>(knot_vector),
                                           std::make_shared<baf::KnotVector>(knot_vector)};
    std::vector<baf::ControlPoint> control_points(8, baf::ControlPoint{2.0});
    control_points[6].SetValue(0, 1.5);
    std::shared_ptr<spl::BSpline<3>> bspline_3d_ptr = std::make_shared<spl::BSpline<3>>(knot_vector_ptr, degree,
                                                                                        control_points);
    bspline_3d_any_ = std::make_any<std::shared_ptr<spl::BSpline<3>>>(bspline_3d_ptr);
  }

 protected:
  std::any bspline_3d_any_;
};

TEST_F(A3DAnyBSplineForAnyCasts, Returns1_5AsSeventhControlPointAfterCastTo3DBSpline) {  // NOLINT
  ASSERT_THAT(util::any_casts::GetSpline<3>(bspline_3d_any_)->GetControlPoint(6).GetValueForDimension(Dimension{0}),
      DoubleEq(1.5));
}

TEST_F(A3DAnyBSplineForAnyCasts, ReturnsSplineDimension3) {  // NOLINT
  ASSERT_THAT(util::any_casts::GetSplineDimension<3>(bspline_3d_any_), 3);
}

TEST_F(A3DAnyBSplineForAnyCasts, ReturnsThatTheSplineIsNotRational) {  // NOLINT
  ASSERT_THAT(util::any_casts::IsRational<3>(bspline_3d_any_), false);
}

class A4DAnyNURBSForAnyCasts : public Test {
 public:
  A4DAnyNURBSForAnyCasts() {
    std::array<Degree, 4> degree = {Degree(1), Degree(1), Degree(1), Degree(1)};
    baf::KnotVector knot_vector({ParametricCoordinate(0), ParametricCoordinate(0), ParametricCoordinate(1),
                                 ParametricCoordinate(1)});
    baf::KnotVectors<4> knot_vector_ptr =
        {std::make_shared<baf::KnotVector>(knot_vector), std::make_shared<baf::KnotVector>(knot_vector),
         std::make_shared<baf::KnotVector>(knot_vector), std::make_shared<baf::KnotVector>(knot_vector)};
    std::vector<baf::ControlPoint> control_points(16, baf::ControlPoint{2.0});
    control_points[15].SetValue(0, 0.4);
    std::vector<double> weights(16, 1.2);
    std::shared_ptr<spl::NURBS<4>> nurbs_4d_ptr = std::make_shared<spl::NURBS<4>>(knot_vector_ptr, degree,
                                                                                  control_points, weights);
    nurbs_4d_any_ = std::make_any<std::shared_ptr<spl::NURBS<4>>>(nurbs_4d_ptr);
  }

 protected:
  std::any nurbs_4d_any_;
};

TEST_F(A4DAnyNURBSForAnyCasts, Returns0_4AsLastControlPointAfterCastTo4DNURBS) {  // NOLINT
  ASSERT_THAT(util::any_casts::GetSpline<4>(nurbs_4d_any_)->GetControlPoint(15).GetValueForDimension(Dimension{0}),
      DoubleEq(0.4));
}

TEST_F(A4DAnyNURBSForAnyCasts, ReturnsSplineDimension4) {  // NOLINT
  ASSERT_THAT(util::any_casts::GetSplineDimension<4>(nurbs_4d_any_), 4);
}

TEST_F(A4DAnyNURBSForAnyCasts, ReturnsThatTheSplineIsRational) {  // NOLINT
  ASSERT_THAT(util::any_casts::IsRational<4>(nurbs_4d_any_), true);
}
