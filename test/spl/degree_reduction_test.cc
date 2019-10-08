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

using namespace splinelib::src;

class BSplineFig5_35ForDegreeElevationAndReductionForDimension0 : public Test {  // NOLINT
 public:
  BSplineFig5_35ForDegreeElevationAndReductionForDimension0() {
    std::array<Degree, 1> degree = {Degree{3}};
    baf::KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(baf::KnotVector(
        {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0},
         ParametricCoordinate{0.3},
         ParametricCoordinate{0.6}, ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{1},
         ParametricCoordinate{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({3.0, 0.0}))
    };
    original_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    elevated_and_reduced_ = std::make_shared<spl::BSpline<1>>(*original_);
    elevated_and_reduced_->ElevateDegreeForDimension(0);
    successful_ = elevated_and_reduced_->ReduceDegreeForDimension(0);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> original_;
  std::shared_ptr<spl::BSpline<1>> elevated_and_reduced_;
  bool successful_;
};

TEST_F(BSplineFig5_35ForDegreeElevationAndReductionForDimension0, // NOLINT
       DegreeReductionWithoutPreviousDegreeElevationFails) { // NOLINT
  std::shared_ptr<spl::BSpline<1>> test_spline = std::make_shared<spl::BSpline<1>>(*original_);
  bool successful = test_spline->ReduceDegreeForDimension(0);
  ASSERT_THAT(successful, false);
  ASSERT_THAT(test_spline->GetDegree(0), original_->GetDegree(0));
  ASSERT_THAT(original_->AreGeometricallyEqual(*elevated_and_reduced_), true);
}

TEST_F(BSplineFig5_35ForDegreeElevationAndReductionForDimension0, DegreeReductionWasSuccessful) {  // NOLINT
  ASSERT_THAT(successful_, true);
}

TEST_F(BSplineFig5_35ForDegreeElevationAndReductionForDimension0, HasUnchangedDegreeInDimension0) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetDegree(0).Get(), elevated_and_reduced_->GetDegree(0).Get());
}

TEST_F(BSplineFig5_35ForDegreeElevationAndReductionForDimension0, HasUnchangedNumberOfKnotsInDimension0) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots());
}

TEST_F(BSplineFig5_35ForDegreeElevationAndReductionForDimension0, HasUnchangedNumberOfControlPoints) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetNumberOfControlPoints(), original_->GetNumberOfControlPoints());
}

TEST_F(BSplineFig5_35ForDegreeElevationAndReductionForDimension0, // NOLINT
       DoesNotChangeGeometricallyAfterDegreeElevationAndReduction) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->AreGeometricallyEqual(*original_), true);
}

class ALinearNURBSForDegreeElevationAndReductionForDimension0 : public Test {  // NOLINT
 public:
  ALinearNURBSForDegreeElevationAndReductionForDimension0() {
    std::array<Degree, 1> degree = {Degree{1}};
    baf::KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(baf::KnotVector(
        {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.3}, ParametricCoordinate{0.6},
         ParametricCoordinate{1},
         ParametricCoordinate{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0}))
    };
    std::vector<double> weights = {1.0, 1.5, 0.5, 1.0};
    original_ = std::make_shared<spl::NURBS<1>>(knot_vector_before, degree, control_points, weights);
    elevated_and_reduced_ = std::make_shared<spl::NURBS<1>>(*original_);
    elevated_and_reduced_->ElevateDegreeForDimension(0);
    successful_ = elevated_and_reduced_->ReduceDegreeForDimension(0);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> original_;
  std::shared_ptr<spl::NURBS<1>> elevated_and_reduced_;
  bool successful_;
};

TEST_F(ALinearNURBSForDegreeElevationAndReductionForDimension0, DegreeReductionWasSuccessful) {  // NOLINT
  ASSERT_THAT(successful_, true);
}

TEST_F(ALinearNURBSForDegreeElevationAndReductionForDimension0, HasUnchangedDegreeInDimension0) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetDegree(0).Get(), original_->GetDegree(0).Get());
}

TEST_F(ALinearNURBSForDegreeElevationAndReductionForDimension0, HasUnchangedNumberOfKnotsInDimension0) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots());
}

TEST_F(ALinearNURBSForDegreeElevationAndReductionForDimension0, HasUnchangedNumberOfControlPoints) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetNumberOfControlPoints(), original_->GetNumberOfControlPoints());
}

TEST_F(ALinearNURBSForDegreeElevationAndReductionForDimension0, // NOLINT
       DoesNotChangeGeometricallyAfterDegreeElevationAndReduction) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->AreGeometricallyEqual(*original_), true);
}

class A2DBSplineForDegreeElevationAndReductionForDimension0 : public Test {  // NOLINT
 public:
  A2DBSplineForDegreeElevationAndReductionForDimension0() {
    std::array<Degree, 2> degree = {Degree{2}, Degree{1}};
    baf::KnotVectors<2> knot_vector_before = {
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.3},
             ParametricCoordinate{0.6},
             ParametricCoordinate{0.6}, ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1}, ParametricCoordinate{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 1.0})), baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})), baf::ControlPoint(std::vector<double>({2.0, 3.0})),
        baf::ControlPoint(std::vector<double>({3.0, 3.0})), baf::ControlPoint(std::vector<double>({3.0, 4.0})),

        baf::ControlPoint(std::vector<double>({2.0, 1.0})), baf::ControlPoint(std::vector<double>({1.5, 2.5})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0})), baf::ControlPoint(std::vector<double>({2.5, 3.5})),
        baf::ControlPoint(std::vector<double>({3.5, 4.0})), baf::ControlPoint(std::vector<double>({5.0, 5.0}))
    };
    original_ = std::make_shared<spl::BSpline<2>>(knot_vector_before, degree, control_points);
    elevated_and_reduced_ = std::make_shared<spl::BSpline<2>>(*original_);
    elevated_and_reduced_->ElevateDegreeForDimension(0);
    successful_ = elevated_and_reduced_->ReduceDegreeForDimension(0);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> original_;
  std::shared_ptr<spl::BSpline<2>> elevated_and_reduced_;
  bool successful_;
};

TEST_F(A2DBSplineForDegreeElevationAndReductionForDimension0, DegreeReductionForDimension0WasSuccessful) {  // NOLINT
  ASSERT_THAT(successful_, true);
}

TEST_F(A2DBSplineForDegreeElevationAndReductionForDimension0, HasUnchangedDegreeInAllDimensions) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetDegree(0).Get(), original_->GetDegree(0).Get());
  ASSERT_THAT(elevated_and_reduced_->GetDegree(1).Get(), original_->GetDegree(1).Get());
}

TEST_F(A2DBSplineForDegreeElevationAndReductionForDimension0, HasUnchangedNumberOfKnotsInAllDimensions) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots());
}

TEST_F(A2DBSplineForDegreeElevationAndReductionForDimension0, // NOLINT
       HasUnchangedNumberOfControlPointsInAllDirections) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0]);
  ASSERT_THAT(elevated_and_reduced_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
}

TEST_F(A2DBSplineForDegreeElevationAndReductionForDimension0, // NOLINT
       DoesNotChangeGeometricallyAfterDegreeElevationAndRedution) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->AreGeometricallyEqual(*original_), true);
}

TEST_F(A2DBSplineForDegreeElevationAndReductionForDimension0, // NOLINT
       DoesNotChangeGeometricallyAfterMoreDegreeElevationsAndReductions) {  // NOLINT
  elevated_and_reduced_->ElevateDegreeForDimension(1);
  elevated_and_reduced_->ElevateDegreeForDimension(0);
  elevated_and_reduced_->ElevateDegreeForDimension(1);
  successful_ = elevated_and_reduced_->ReduceDegreeForDimension(0);
  successful_ = successful_ && elevated_and_reduced_->ReduceDegreeForDimension(1);
  successful_ = successful_ && elevated_and_reduced_->ReduceDegreeForDimension(1);
  ASSERT_THAT(successful_, true);
  ASSERT_THAT(elevated_and_reduced_->AreGeometricallyEqual(*original_), true);
}

class Random2DNURBSForDegreeElevationAndReductionForDimension0 : public Test {  // NOLINT
 public:
  Random2DNURBSForDegreeElevationAndReductionForDimension0() {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0.0}, ParametricCoordinate{1.0}};
    spl::RandomNURBSGenerator<2> spline_generator(limits, 3, 3);
    original_ = std::make_shared<spl::NURBS<2>>(spline_generator);
    elevated_and_reduced_ = std::make_shared<spl::NURBS<2>>(*original_);
    elevated_and_reduced_->ElevateDegreeForDimension(0);
    successful_ = elevated_and_reduced_->ReduceDegreeForDimension(0);
  }

 protected:
  std::shared_ptr<spl::NURBS<2>> original_;
  std::shared_ptr<spl::NURBS<2>> elevated_and_reduced_;
  bool successful_;
};

TEST_F(Random2DNURBSForDegreeElevationAndReductionForDimension0, DegreeReductionForDimension0WasSuccessful) {  // NOLINT
  ASSERT_THAT(successful_, true);
}

TEST_F(Random2DNURBSForDegreeElevationAndReductionForDimension0, HasUnchangedDegreeInAllDimensions) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetDegree(0).Get(), original_->GetDegree(0).Get());
  ASSERT_THAT(elevated_and_reduced_->GetDegree(1).Get(), original_->GetDegree(1).Get());
}

TEST_F(Random2DNURBSForDegreeElevationAndReductionForDimension0, HasUnchangedNumberOfKnotsInAllDimensions) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots());
}

TEST_F(Random2DNURBSForDegreeElevationAndReductionForDimension0, // NOLINT
       HasUnchangedNumberOfControlPointsInAllDimensions) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0]);
  ASSERT_THAT(elevated_and_reduced_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
}

TEST_F(Random2DNURBSForDegreeElevationAndReductionForDimension0, // NOLINT
       DoesNotChangeGeometricallyAfterDegreeElevationAndReduction) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->AreGeometricallyEqual(*original_), true);
}

class Random3DBSplineForDegreeElevationAndReductionForDimension0 : public Test {  // NOLINT
 public:
  Random3DBSplineForDegreeElevationAndReductionForDimension0() {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0.0}, ParametricCoordinate{1.0}};
    spl::RandomBSplineGenerator<3> spline_generator(limits, 4, 3);
    original_ = std::make_shared<spl::BSpline<3>>(spline_generator);
    elevated_and_reduced_ = std::make_shared<spl::BSpline<3>>(*original_);
    elevated_and_reduced_->ElevateDegreeForDimension(0);
    successful_ = elevated_and_reduced_->ReduceDegreeForDimension(0);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> original_;
  std::shared_ptr<spl::BSpline<3>> elevated_and_reduced_;
  bool successful_;
};

TEST_F(Random3DBSplineForDegreeElevationAndReductionForDimension0, // NOLINT
       DegreeReductionForDimension0WasSuccessful) {  // NOLINT
  ASSERT_THAT(successful_, true);
}

TEST_F(Random3DBSplineForDegreeElevationAndReductionForDimension0, HasUnchangedDegreeInAllDimensions) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetDegree(0).Get(), original_->GetDegree(0).Get());
  ASSERT_THAT(elevated_and_reduced_->GetDegree(1).Get(), original_->GetDegree(1).Get());
  ASSERT_THAT(elevated_and_reduced_->GetDegree(2).Get(), original_->GetDegree(2).Get());
}

TEST_F(Random3DBSplineForDegreeElevationAndReductionForDimension0, // NOLINT
       HasUnchangedNumberOfKnotsInAllDimensions) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots());
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(2)->GetNumberOfKnots(),
              original_->GetKnotVector(2)->GetNumberOfKnots());
}

TEST_F(Random3DBSplineForDegreeElevationAndReductionForDimension0, // NOLINT
       HasUnchangedNumberOfControlPointsInAllDimensions) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0]);
  ASSERT_THAT(elevated_and_reduced_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
}

TEST_F(Random3DBSplineForDegreeElevationAndReductionForDimension0, // NOLINT
       DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->AreGeometricallyEqual(*original_), true);
}

TEST_F(Random3DBSplineForDegreeElevationAndReductionForDimension0, // NOLINT
       DoesNotChangeGeometricallyAfterMoreDegreeElevationsAndReductions) {  // NOLINT
  elevated_and_reduced_->ElevateDegreeForDimension(2);
  elevated_and_reduced_->ElevateDegreeForDimension(2);
  elevated_and_reduced_->ElevateDegreeForDimension(1);
  elevated_and_reduced_->ElevateDegreeForDimension(0);
  successful_ = elevated_and_reduced_->ReduceDegreeForDimension(1);
  successful_ = successful_ && elevated_and_reduced_->ReduceDegreeForDimension(0);
  successful_ = successful_ && elevated_and_reduced_->ReduceDegreeForDimension(2);
  successful_ = successful_ && elevated_and_reduced_->ReduceDegreeForDimension(2);
  ASSERT_THAT(successful_, true);
  ASSERT_THAT(elevated_and_reduced_->AreGeometricallyEqual(*original_), true);
}

class Random3DNURBSForDegreeElevationAndReductionForDimension1 : public Test {  // NOLINT
 public:
  Random3DNURBSForDegreeElevationAndReductionForDimension1() {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0.0}, ParametricCoordinate{1.0}};
    spl::RandomNURBSGenerator<3> spline_generator(limits, 3, 3);
    original_ = std::make_shared<spl::NURBS<3>>(spline_generator);
    elevated_and_reduced_ = std::make_shared<spl::NURBS<3>>(*original_);
    elevated_and_reduced_->ElevateDegreeForDimension(1);
    successful_ = elevated_and_reduced_->ReduceDegreeForDimension(1);
  }

 protected:
  std::shared_ptr<spl::NURBS<3>> original_;
  std::shared_ptr<spl::NURBS<3>> elevated_and_reduced_;
  bool successful_;
};

TEST_F(Random3DNURBSForDegreeElevationAndReductionForDimension1, DegreeReductionForDimension1WasSuccessful) {  // NOLINT
  ASSERT_THAT(successful_, true);
}

TEST_F(Random3DNURBSForDegreeElevationAndReductionForDimension1, HasUnchangedDegreeInAllDimensions) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetDegree(0).Get(), original_->GetDegree(0).Get());
  ASSERT_THAT(elevated_and_reduced_->GetDegree(1).Get(), original_->GetDegree(1).Get());
  ASSERT_THAT(elevated_and_reduced_->GetDegree(2).Get(), original_->GetDegree(2).Get());
}

TEST_F(Random3DNURBSForDegreeElevationAndReductionForDimension1, // NOLINT
       HasUnchangedNumberOfKnotsInAllDimensions) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(1)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(1)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots());
  ASSERT_THAT(elevated_and_reduced_->GetKnotVector(2)->GetNumberOfKnots(),
              original_->GetKnotVector(2)->GetNumberOfKnots());
}

TEST_F(Random3DNURBSForDegreeElevationAndReductionForDimension1, // NOLINT
       HasUnchangedNumberOfControlPointsInAllDimensions) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0]);
  ASSERT_THAT(elevated_and_reduced_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
  ASSERT_THAT(elevated_and_reduced_->GetPointsPerDirection()[2], original_->GetPointsPerDirection()[2]);
}

TEST_F(Random3DNURBSForDegreeElevationAndReductionForDimension1, // NOLINT
       DoesNotChangeGeometricallyAfterDegreeElevationAndReduction) {  // NOLINT
  ASSERT_THAT(elevated_and_reduced_->AreGeometricallyEqual(*original_), true);
}

TEST_F(Random3DNURBSForDegreeElevationAndReductionForDimension1, // NOLINT
       DoesNotChangeGeometricallyAfterMoreDegreeElevationsAndReductions) {  // NOLINT
  elevated_and_reduced_->ElevateDegreeForDimension(2);
  elevated_and_reduced_->ElevateDegreeForDimension(2);
  elevated_and_reduced_->ElevateDegreeForDimension(1);
  elevated_and_reduced_->ElevateDegreeForDimension(0);
  successful_ = elevated_and_reduced_->ReduceDegreeForDimension(1);
  successful_ = successful_ && elevated_and_reduced_->ReduceDegreeForDimension(0);
  successful_ = successful_ && elevated_and_reduced_->ReduceDegreeForDimension(2);
  successful_ = successful_ && elevated_and_reduced_->ReduceDegreeForDimension(2);
  ASSERT_THAT(successful_, true);
  ASSERT_THAT(elevated_and_reduced_->AreGeometricallyEqual(*original_), true);
}
