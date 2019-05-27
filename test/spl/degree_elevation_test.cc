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
#include "random_b_spline_generator.h"
#include "random_nurbs_generator.h"

using testing::Test;

template<int DIM>
void PrintSpline(std::shared_ptr<spl::BSpline<DIM>> spline) {
  std::cout << std::endl << "degrees:" << std::endl;
  for (int i = 0; i < DIM; ++i) {
    std::cout << spline->GetDegree(i).get() << "   ";
  }
  std::cout << std::endl << std::endl << "knot vectors:" << std::endl;
  for (int i = 0; i < DIM; ++i) {
    auto kv = spline->GetKnotVector(i);
    for (const auto &knot : *kv) {
      std::cout << knot.get() << "  ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl << "control points:" << std::endl;
  for (int i = 0; i < spline->GetNumberOfControlPoints(); ++i) {
    for (int j = 0; j < spline->GetPointDim(); ++j) {
      std::cout << spline->GetControlPoint({i}, j) << "  ";
    }
    std::cout << std::endl;
  }
}

class BSpline1DFig5_35 : public Test {  // NOLINT
 public:
  BSpline1DFig5_35() {
    std::array<Degree, 1> degree = {Degree{3}};
    KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6},
                         ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    KnotVectors<1> knot_vector_after = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6},
                         ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({3.0, 0.0}))
    };
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(knot_vector_after, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(BSpline1DFig5_35, ElevatesDegreeFrom3To4Correctly) {  // NOLINT
  bspline_1d_after_->ElevateDegree(0);
  ASSERT_THAT(bspline_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 4);
  ASSERT_THAT(bspline_1d_after_->GetNumberOfControlPoints(), bspline_1d_before_->GetNumberOfControlPoints() + 3);
  ASSERT_THAT(bspline_1d_after_->AreGeometricallyEqual(*bspline_1d_before_), true);
}

class ALinearBSpline : public Test {  // NOLINT
 public:
  ALinearBSpline() {
    std::array<Degree, 1> degree = {Degree{1}};
    KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{1},
                         ParamCoord{1}}))};
    KnotVectors<1> knot_vector_after = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{1},
                         ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0}))
    };
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(knot_vector_after, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(ALinearBSpline, ElevatesDegreeFrom1To2Correctly) {  // NOLINT
  std::array<Degree, 1> degree = {Degree{2}};
  KnotVectors<1> knot_vector = {std::make_shared<baf::KnotVector>(
      baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.3}, ParamCoord{0.6},
                       ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
  std::vector<baf::ControlPoint> control_points = {
      baf::ControlPoint(std::vector<double>({1.0, 1.0})),
      baf::ControlPoint(std::vector<double>({1.0, 1.5})),
      baf::ControlPoint(std::vector<double>({1.0, 2.0})),
      baf::ControlPoint(std::vector<double>({1.5, 2.0})),
      baf::ControlPoint(std::vector<double>({2.0, 2.0})),
      baf::ControlPoint(std::vector<double>({2.0, 2.5})),
      baf::ControlPoint(std::vector<double>({2.0, 3.0}))
  };
  auto test = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);

  bspline_1d_after_->ElevateDegree(0);
  ASSERT_THAT(bspline_1d_after_->AreEqual(*test), true);
  ASSERT_THAT(bspline_1d_before_->AreGeometricallyEqual(*test), true);
  ASSERT_THAT(bspline_1d_after_->AreGeometricallyEqual(*test), true);
  ASSERT_THAT(bspline_1d_after_->AreGeometricallyEqual(*bspline_1d_before_), true);
}

class AQuadraticBSpline : public Test {  // NOLINT
 public:
  AQuadraticBSpline() {
    std::array<Degree, 1> degree = {Degree{2}};
    KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{0.6},
                         ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    KnotVectors<1> knot_vector_after = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{0.6},
                         ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0})),
        baf::ControlPoint(std::vector<double>({3.0, 3.0})),
        baf::ControlPoint(std::vector<double>({3.0, 4.0}))
    };
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(knot_vector_after, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(AQuadraticBSpline, ElevatesDegreeFrom2To3Correctly) {  // NOLINT
  std::array<Degree, 1> degree = {Degree{3}};
  KnotVectors<1> knot_vector = {std::make_shared<baf::KnotVector>(
      baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.3},
                       ParamCoord{0.6}, ParamCoord{0.6}, ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1},
                       ParamCoord{1}}))};
  std::vector<baf::ControlPoint> control_points = {
      baf::ControlPoint(std::vector<double>({1.0, 1.0})),
      baf::ControlPoint(std::vector<double>({1.0, 5.0 / 3.0})),
      baf::ControlPoint(std::vector<double>({7.0 / 6.0, 2.0})),
      baf::ControlPoint(std::vector<double>({11.0 / 6.0, 2.0})),
      baf::ControlPoint(std::vector<double>({2.0, 7.0 / 3.0})),
      baf::ControlPoint(std::vector<double>({2.0, 3.0})),
      baf::ControlPoint(std::vector<double>({8.0 / 3.0, 3.0})),
      baf::ControlPoint(std::vector<double>({3.0, 10.0 / 3.0})),
      baf::ControlPoint(std::vector<double>({3.0, 4.0}))
  };
  auto test = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);

  bspline_1d_after_->ElevateDegree(0);
  ASSERT_THAT(bspline_1d_after_->AreEqual(*test), true);
  ASSERT_THAT(bspline_1d_before_->AreGeometricallyEqual(*test), true);
  ASSERT_THAT(bspline_1d_after_->AreGeometricallyEqual(*test), true);
  ASSERT_THAT(bspline_1d_after_->AreGeometricallyEqual(*bspline_1d_before_), true);
}

class Random1DBSplineForDegreeElevation : public Test {  // NOLINT
 public:
  Random1DBSplineForDegreeElevation() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomBSplineGenerator<1> spline_generator(limits, 10, 3);
    spl::BSpline<1> b_spline(spline_generator);
    original_ = std::make_shared<spl::BSpline<1>>(b_spline);

    spl::BSpline<1> elevation_spline(b_spline);
    elevation_spline.ElevateDegree(0);
    after_elevation = std::make_shared<spl::BSpline<1>>(elevation_spline);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> original_;
  std::shared_ptr<spl::BSpline<1>> after_elevation;
};

TEST_F(Random1DBSplineForDegreeElevation, HasELevatedDegree) {
  ASSERT_THAT(after_elevation->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
}

TEST_F(Random1DBSplineForDegreeElevation, HasMoreKnots) {
  ASSERT_THAT(after_elevation->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots()
                  + original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
}

TEST_F(Random1DBSplineForDegreeElevation, HasMoreControlPoints) {
  ASSERT_THAT(after_elevation->GetNumberOfControlPoints(),
              original_->GetNumberOfControlPoints() + original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
}

TEST_F(Random1DBSplineForDegreeElevation, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(after_elevation->AreGeometricallyEqual(*original_), true);
}
