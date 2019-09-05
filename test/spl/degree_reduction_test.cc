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
#include "random_b_spline_generator.h"
#include "nurbs.h"
#include "vtk_writer.h"

using testing::Test;

class ALinearBSpline : public Test {  // NOLINT
 public:
  ALinearBSpline() {
    std::array<Degree, 1> degree = {Degree{1}};
    KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{1},
                         ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0}))
    };
    b_spline_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    b_spline_after_ = std::make_shared<spl::BSpline<1>>(*b_spline_before_);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> b_spline_before_;
  std::shared_ptr<spl::BSpline<1>> b_spline_after_;
};

TEST_F(ALinearBSpline, DoesNotChangeGeometricallyAfterDegreeReduction) {  // NOLINT
  b_spline_after_->ElevateDegreeForDimension(0);
  b_spline_after_->ElevateDegreeForDimension(0);
  bool successful = b_spline_after_->ReduceDegreeForDimension(0, 0.0001);
  successful = successful && b_spline_after_->ReduceDegreeForDimension(0, 0.0001);
  ASSERT_THAT(successful, true);
  ASSERT_THAT(b_spline_after_->AreGeometricallyEqual(*b_spline_before_), true);
}

class BSpline1DFig5_39 : public Test {  // NOLINT
 public:
  BSpline1DFig5_39() {
    std::array<Degree, 1> degree = {Degree{4}};
    KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3},
                         ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1},
                         ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({-0.5, 0.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 0.5})),
        baf::ControlPoint(std::vector<double>({-2.0, 1.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 1.5})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.5})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.5})),
        baf::ControlPoint(std::vector<double>({0.5, 0.0}))
    };
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(BSpline1DFig5_39, ReducesDegreeCorrectly) {  // NOLINT
  bspline_1d_after_->ElevateDegreeForDimension(0);
  bool successful = bspline_1d_after_->ReduceDegreeForDimension(0, 0.0001);
  ASSERT_THAT(successful, true);
  ASSERT_THAT(bspline_1d_after_->GetDegree(0).get(), bspline_1d_before_->GetDegree(0).get());
  ASSERT_THAT(bspline_1d_after_->AreGeometricallyEqual(*bspline_1d_before_), true);
}

class A2DBSplineForDegreeReduction : public Test {  // NOLINT
 public:
  A2DBSplineForDegreeReduction() {
    std::array<Degree, 2> degree = {Degree{2}, Degree{1}};
    KnotVectors<2> knot_vector_before = {
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6},
                             ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})), baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})), baf::ControlPoint(std::vector<double>({0.0, 3.0})),
        baf::ControlPoint(std::vector<double>({0.0, 4.0})), baf::ControlPoint(std::vector<double>({0.0, 5.0})),

        baf::ControlPoint(std::vector<double>({2.0, 0.0})), baf::ControlPoint(std::vector<double>({2.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})), baf::ControlPoint(std::vector<double>({2.0, 3.0})),
        baf::ControlPoint(std::vector<double>({2.0, 4.0})), baf::ControlPoint(std::vector<double>({2.0, 5.0}))
    };
    original_ = std::make_shared<spl::BSpline<2>>(knot_vector_before, degree, control_points);
    after_reduction_ = std::make_shared<spl::BSpline<2>>(*original_);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> original_;
  std::shared_ptr<spl::BSpline<2>> after_reduction_;
};

TEST_F(A2DBSplineForDegreeReduction, ReducesDegreeCorrectly) {  // NOLINT
  after_reduction_->ElevateDegreeForDimension(0);
  bool successful = after_reduction_->ReduceDegreeForDimension(0, 0.001);
  ASSERT_THAT(successful, true);
  ASSERT_THAT(after_reduction_->GetDegree(1).get(), original_->GetDegree(1).get());
  ASSERT_THAT(after_reduction_->AreGeometricallyEqual(*original_), true);
}

class Random3DBSplineForDegreeReduction : public Test {  // NOLINT
 public:
  Random3DBSplineForDegreeReduction() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomBSplineGenerator<3> spline_generator(limits, 4, 3);
    spl::BSpline<3> b_spline(spline_generator);
    b_spline_before_ = std::make_shared<spl::BSpline<3>>(b_spline);
    b_spline_after_ = std::make_shared<spl::BSpline<3>>(b_spline);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> b_spline_before_;
  std::shared_ptr<spl::BSpline<3>> b_spline_after_;
};

TEST_F(Random3DBSplineForDegreeReduction, ReducesDegreeCorrectly) {  // NOLINT
  b_spline_after_->ElevateDegreeForDimension(0);
  b_spline_after_->ElevateDegreeForDimension(1);
  b_spline_after_->ElevateDegreeForDimension(2);
  b_spline_after_->ReduceDegreeForDimension(0, 1e-10);
  b_spline_after_->ReduceDegreeForDimension(1, 1e-10);
  b_spline_after_->ReduceDegreeForDimension(2, 1e-10);
  ASSERT_THAT(b_spline_after_->AreGeometricallyEqual(*b_spline_before_), true);
}
