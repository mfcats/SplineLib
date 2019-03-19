/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "any_casts.h"

#include "gmock/gmock.h"

using testing::Test;
using testing::DoubleEq;

class AnySplines : public Test {
 public:
  AnySplines() {
    std::array<Degree, 1> degree = {Degree(1)};
    KnotVectors<1> knot_vector_ptr =
        {std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord(0), ParamCoord(0), ParamCoord(1),
                                                            ParamCoord(1)}))};
    std::vector<baf::ControlPoint> control_points = {baf::ControlPoint{0.0}, baf::ControlPoint{1.2}};
    spl::BSpline<1> b_spline_1d(knot_vector_ptr, degree, control_points);
    std::shared_ptr<spl::BSpline<1>> b_spline_1d_ptr = std::make_shared<spl::BSpline<1>>(b_spline_1d);
    b_spline_1d_any_ = std::make_any<std::shared_ptr<spl::BSpline<1>>>(b_spline_1d_ptr);

    std::array<Degree, 2> degrees = {Degree(1), Degree(1)};
    KnotVectors<2> knot_vector_ptrs = {knot_vector_ptr[0], knot_vector_ptr[0]};
    control_points = {baf::ControlPoint{0.0}, baf::ControlPoint{1.0}, baf::ControlPoint{2.0}, baf::ControlPoint{3.0}};
    std::vector<double> weights = {0.9, 0.8, 1.0, 1.2};
    spl::NURBS<2> nurbs_2d(knot_vector_ptrs, degrees, control_points, weights);
    std::shared_ptr<spl::NURBS<2>> nurbs_2d_ptr = std::make_shared<spl::NURBS<2>>(nurbs_2d);
    nurbs_2d_any_ = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbs_2d_ptr);
  }

  std::any b_spline_1d_any_;
  std::any nurbs_2d_any_;
};

TEST_F(AnySplines, CanBeCheckedForSplineDimension) {  // NOLINT
  ASSERT_THAT(util::AnyCasts::GetSplineDimension(b_spline_1d_any_), 1);
  ASSERT_THAT(util::AnyCasts::GetSplineDimension(nurbs_2d_any_), 2);
  ASSERT_THROW(util::AnyCasts::GetSplineDimension(std::make_any<int>(8)), std::runtime_error);
}

TEST_F(AnySplines, CanBeCastedToSplines) {  // NOLINT
  ASSERT_THAT(util::AnyCasts::GetSpline<1>(b_spline_1d_any_)->GetControlPoint({1}).GetValue(0), DoubleEq(1.2));
  ASSERT_THAT(util::AnyCasts::GetSpline<2>(nurbs_2d_any_)->GetWeights()[1], DoubleEq(0.8));
  ASSERT_THROW(util::AnyCasts::GetSpline<1>(std::make_any<int>(8)), std::runtime_error);
}

TEST_F(AnySplines, CanBeCheckedIfRational) {  // NOLINT
  ASSERT_THAT(util::AnyCasts::IsRational<1>(b_spline_1d_any_), false);
  ASSERT_THAT(util::AnyCasts::IsRational<2>(nurbs_2d_any_), true);
  ASSERT_THROW(util::AnyCasts::IsRational<1>(std::make_any<int>(8)), std::runtime_error);
}
