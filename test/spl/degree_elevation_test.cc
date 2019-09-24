/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <memory>

#include "gmock/gmock.h"

#include "b_spline.h"
#include "nurbs.h"
#include "random_b_spline_generator.h"
#include "random_nurbs_generator.h"

using testing::Test;

class BSplineFig5_35ForDegreeElevationForDimension0 : public Test {  // NOLINT
 public:
  BSplineFig5_35ForDegreeElevationForDimension0() {
    std::array<Degree, 1> degree = {Degree{3}};
    KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(
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
    original_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    elevated_ = std::make_shared<spl::BSpline<1>>(*original_);
    elevated_->ElevateDegreeForDimension(0);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> original_;
  std::shared_ptr<spl::BSpline<1>> elevated_;
};

TEST_F(BSplineFig5_35ForDegreeElevationForDimension0, HasELevatedDegreeInDimension0) {  // NOLINT
  ASSERT_THAT(elevated_->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
}

TEST_F(BSplineFig5_35ForDegreeElevationForDimension0, HasNumberOfDifferentKnotsMoreKnots) {  // NOLINT
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots() +
                  original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
}

TEST_F(BSplineFig5_35ForDegreeElevationForDimension0, HasNumberOfNonZeroKnotSpansMoreControlPoints) {  // NOLINT
  ASSERT_THAT(elevated_->GetNumberOfControlPoints(), original_->GetNumberOfControlPoints() +
      original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
}

TEST_F(BSplineFig5_35ForDegreeElevationForDimension0, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
}

class ALinearNURBSForDegreeElevationForDimension0 : public Test {  // NOLINT
 public:
  ALinearNURBSForDegreeElevationForDimension0() {
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
    std::vector<double> weights = {1.0, 1.5, 0.5, 1.0};
    original_ = std::make_shared<spl::NURBS<1>>(knot_vector_before, degree, control_points, weights);
    elevated_ = std::make_shared<spl::NURBS<1>>(*original_);
    elevated_->ElevateDegreeForDimension(0);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> original_;
  std::shared_ptr<spl::NURBS<1>> elevated_;
};

TEST_F(ALinearNURBSForDegreeElevationForDimension0, HasELevatedDegreeInDimension0) {  // NOLINT
  ASSERT_THAT(elevated_->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
}

TEST_F(ALinearNURBSForDegreeElevationForDimension0, HasNumberOfDifferentKnotsMoreKnots) {  // NOLINT
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots() +
                  original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
}

TEST_F(ALinearNURBSForDegreeElevationForDimension0, HasNumberOfNonZeroKnotSpansMoreControlPoints) {  // NOLINT
  ASSERT_THAT(elevated_->GetNumberOfControlPoints(), original_->GetNumberOfControlPoints() +
      original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
}

TEST_F(ALinearNURBSForDegreeElevationForDimension0, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
}

class A2DBSplineForDegreeElevationForDimension0 : public Test {  // NOLINT
 public:
  A2DBSplineForDegreeElevationForDimension0() {
    std::array<Degree, 2> degree = {Degree{2}, Degree{1}};
    KnotVectors<2> knot_vector_before = {
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6},
                             ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 1.0})), baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})), baf::ControlPoint(std::vector<double>({2.0, 3.0})),
        baf::ControlPoint(std::vector<double>({3.0, 3.0})), baf::ControlPoint(std::vector<double>({3.0, 4.0})),

        baf::ControlPoint(std::vector<double>({2.0, 1.0})), baf::ControlPoint(std::vector<double>({1.5, 2.5})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0})), baf::ControlPoint(std::vector<double>({2.5, 3.5})),
        baf::ControlPoint(std::vector<double>({3.5, 4.0})), baf::ControlPoint(std::vector<double>({5.0, 5.0}))
    };
    original_ = std::make_shared<spl::BSpline<2>>(knot_vector_before, degree, control_points);
    elevated_ = std::make_shared<spl::BSpline<2>>(*original_);
    elevated_->ElevateDegreeForDimension(0);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> original_;
  std::shared_ptr<spl::BSpline<2>> elevated_;
};

TEST_F(A2DBSplineForDegreeElevationForDimension0, HasELevatedDegreeInDimension0) {  // NOLINT
  ASSERT_THAT(elevated_->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
  ASSERT_THAT(elevated_->GetDegree(1).get(), original_->GetDegree(1).get());
}

TEST_F(A2DBSplineForDegreeElevationForDimension0, HasNumberOfDifferentKnotsMoreKnots) {  // NOLINT
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots() +
                  original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_->GetKnotVector(1)->GetNumberOfKnots(), original_->GetKnotVector(1)->GetNumberOfKnots());
}

TEST_F(A2DBSplineForDegreeElevationForDimension0, HasNumberOfNonZeroKnotSpansMoreControlPoints) {  // NOLINT
  ASSERT_THAT(elevated_->GetPointsPerDirection()[0],
              original_->GetPointsPerDirection()[0] + original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
  ASSERT_THAT(elevated_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
}

TEST_F(A2DBSplineForDegreeElevationForDimension0, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
}

TEST_F(A2DBSplineForDegreeElevationForDimension0, DoesNotChangeGeometricallyAfterMoreDegreeElevations) {  // NOLINT
  elevated_->ElevateDegreeForDimension(1);
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
  elevated_->ElevateDegreeForDimension(0);
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
  elevated_->ElevateDegreeForDimension(1);
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_, 1e-10), true);
}
class Random2DNURBSForDegreeElevationForDimension0 : public Test {  // NOLINT
 public:
  Random2DNURBSForDegreeElevationForDimension0() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomNURBSGenerator<2> spline_generator(limits, 3, 3);
    spl::NURBS<2> b_spline(spline_generator);
    original_ = std::make_shared<spl::NURBS<2>>(b_spline);

    spl::NURBS<2> elevation_spline(b_spline);
    elevation_spline.ElevateDegreeForDimension(0);
    elevated_ = std::make_shared<spl::NURBS<2>>(elevation_spline);
  }

 protected:
  std::shared_ptr<spl::NURBS<2>> original_;
  std::shared_ptr<spl::NURBS<2>> elevated_;
};

TEST_F(Random2DNURBSForDegreeElevationForDimension0, HasELevatedDegreeInDimension0) {  // NOLINT
  ASSERT_THAT(elevated_->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
  ASSERT_THAT(elevated_->GetDegree(1).get(), original_->GetDegree(1).get());
}

TEST_F(Random2DNURBSForDegreeElevationForDimension0, HasNumberOfDifferentKnotsMoreKnots) {  // NOLINT
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots() +
                  original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_->GetKnotVector(1)->GetNumberOfKnots(), original_->GetKnotVector(1)->GetNumberOfKnots());
}

TEST_F(Random2DNURBSForDegreeElevationForDimension0, HasNumberOfNonZeroKnotSpansMoreControlPoints) {  // NOLINT
  ASSERT_THAT(elevated_->GetPointsPerDirection()[0],
              original_->GetPointsPerDirection()[0] + original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
  ASSERT_THAT(elevated_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
}

TEST_F(Random2DNURBSForDegreeElevationForDimension0, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
}

class Random3DBSplineForDegreeElevationForDimension0 : public Test {  // NOLINT
 public:
  Random3DBSplineForDegreeElevationForDimension0() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomBSplineGenerator<3> spline_generator(limits, 4, 3);
    spl::BSpline<3> b_spline(spline_generator);
    original_ = std::make_shared<spl::BSpline<3>>(b_spline);

    spl::BSpline<3> elevation_spline(b_spline);
    elevation_spline.ElevateDegreeForDimension(0);
    elevated_ = std::make_shared<spl::BSpline<3>>(elevation_spline);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> original_;
  std::shared_ptr<spl::BSpline<3>> elevated_;
};

TEST_F(Random3DBSplineForDegreeElevationForDimension0, HasELevatedDegreeInDimension0) {  // NOLINT
  ASSERT_THAT(elevated_->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
  ASSERT_THAT(elevated_->GetDegree(1).get(), original_->GetDegree(1).get());
  ASSERT_THAT(elevated_->GetDegree(2).get(), original_->GetDegree(2).get());
}

TEST_F(Random3DBSplineForDegreeElevationForDimension0, HasNumberOfDifferentKnotsMoreKnots) {  // NOLINT
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots() +
                  original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_->GetKnotVector(1)->GetNumberOfKnots(), original_->GetKnotVector(1)->GetNumberOfKnots());
  ASSERT_THAT(elevated_->GetKnotVector(2)->GetNumberOfKnots(), original_->GetKnotVector(2)->GetNumberOfKnots());
}

TEST_F(Random3DBSplineForDegreeElevationForDimension0, HasNumberOfNonZeroKnotSpansMoreControlPoints) {  // NOLINT
  ASSERT_THAT(elevated_->GetPointsPerDirection()[0],
              original_->GetPointsPerDirection()[0] + original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
  ASSERT_THAT(elevated_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
}

TEST_F(Random3DBSplineForDegreeElevationForDimension0, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
}

TEST_F(Random3DBSplineForDegreeElevationForDimension0, DoesNotChangeGeometricallyAfterMoreDegreeElevations) {  // NOLINT
  elevated_->ElevateDegreeForDimension(1);
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
  elevated_->ElevateDegreeForDimension(2);
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
}

class Random3DNURBSForDegreeElevationForDimension1 : public Test {  // NOLINT
 public:
  Random3DNURBSForDegreeElevationForDimension1() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomNURBSGenerator<3> spline_generator(limits, 3, 3);
    spl::NURBS<3> b_spline(spline_generator);
    original_ = std::make_shared<spl::NURBS<3>>(b_spline);

    spl::NURBS<3> elevation_spline(b_spline);
    elevation_spline.ElevateDegreeForDimension(1);
    elevated_ = std::make_shared<spl::NURBS<3>>(elevation_spline);
  }

 protected:
  std::shared_ptr<spl::NURBS<3>> original_;
  std::shared_ptr<spl::NURBS<3>> elevated_;
};

TEST_F(Random3DNURBSForDegreeElevationForDimension1, HasELevatedDegreeInDimension1) {  // NOLINT
  ASSERT_THAT(elevated_->GetDegree(0).get(), original_->GetDegree(0).get());
  ASSERT_THAT(elevated_->GetDegree(1).get(), original_->GetDegree(1).get() + 1);
  ASSERT_THAT(elevated_->GetDegree(2).get(), original_->GetDegree(2).get());
}

TEST_F(Random3DNURBSForDegreeElevationForDimension1, HasNumberOfDifferentKnotsMoreKnots) {  // NOLINT
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfKnots(), original_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(elevated_->GetKnotVector(1)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(1)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots() +
                  original_->GetKnotVector(1)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_->GetKnotVector(2)->GetNumberOfKnots(), original_->GetKnotVector(2)->GetNumberOfKnots());
}

TEST_F(Random3DNURBSForDegreeElevationForDimension1, HasNumberOfNonZeroKnotSpansMoreControlPoints) {  // NOLINT
  ASSERT_THAT(elevated_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0]);
  ASSERT_THAT(elevated_->GetPointsPerDirection()[1],
              original_->GetPointsPerDirection()[1] + original_->GetKnotVector(1)->GetNumberOfDifferentKnots() - 1);
}

TEST_F(Random3DNURBSForDegreeElevationForDimension1, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
}
