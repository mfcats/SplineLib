/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include <memory>

#include "gmock/gmock.h"

#include "src/spl/b_spline.h"
#include "src/spl/nurbs.h"
#include "test/spl/random/random_spline_utils.h"
#include "test/spl/random/write_random_spline.h"

using testing::Test;

using namespace splinelib::src;

class BSplineFig5_35ForDegreeElevationForDimension0 : public Test {  // NOLINT
 public:
  BSplineFig5_35ForDegreeElevationForDimension0() {
    std::array<Degree, 1> degree = {Degree{3}};
    baf::KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(baf::KnotVector(
        {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0},
         ParametricCoordinate{0.3},
         ParametricCoordinate{0.6}, ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{1},
         ParametricCoordinate{1}}))};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint(std::vector<double>({1.0, 0.0})),
        spl::ControlPoint(std::vector<double>({0.0, 2.0})),
        spl::ControlPoint(std::vector<double>({1.0, 2.0})),
        spl::ControlPoint(std::vector<double>({3.0, 2.0})),
        spl::ControlPoint(std::vector<double>({4.0, 1.0})),
        spl::ControlPoint(std::vector<double>({3.0, 0.0}))
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
  ASSERT_THAT(elevated_->GetDegree(0).Get(), original_->GetDegree(0).Get() + 1);
}

TEST_F(BSplineFig5_35ForDegreeElevationForDimension0, HasNumberOfDifferentKnotsMoreKnots) {  // NOLINT
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots() +
                  original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
}

TEST_F(BSplineFig5_35ForDegreeElevationForDimension0, HasNumberOfNonZeroKnotSpansMoreControlPoints) {  // NOLINT
  ASSERT_THAT(elevated_->GetTotalNumberOfControlPoints(), original_->GetTotalNumberOfControlPoints() +
      original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
}

TEST_F(BSplineFig5_35ForDegreeElevationForDimension0, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
}

TEST_F(BSplineFig5_35ForDegreeElevationForDimension0, WriteRandomSplineToXMLIfAnyPreviousTestFailed) {  // NOLINT
  if (testing::UnitTest::GetInstance()->current_test_case()->failed_test_count() > 0)
    splinelib::test::random_spline_writer::WriteToXML<1>({original_, elevated_}, testing::UnitTest::GetInstance());
}

class ALinearNURBSForDegreeElevationForDimension0 : public Test {  // NOLINT
 public:
  ALinearNURBSForDegreeElevationForDimension0() {
    std::array<Degree, 1> degree = {Degree{1}};
    baf::KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(baf::KnotVector(
        {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.3}, ParametricCoordinate{0.6},
         ParametricCoordinate{1},
         ParametricCoordinate{1}}))};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint(std::vector<double>({1.0, 1.0})),
        spl::ControlPoint(std::vector<double>({1.0, 2.0})),
        spl::ControlPoint(std::vector<double>({2.0, 2.0})),
        spl::ControlPoint(std::vector<double>({2.0, 3.0}))
    };
    std::vector<Weight> weights = {Weight{1.0}, Weight{1.5}, Weight{0.5}, Weight{1.0}};
    original_ = std::make_shared<spl::NURBS<1>>(knot_vector_before, degree, control_points, weights);
    elevated_ = std::make_shared<spl::NURBS<1>>(*original_);
    elevated_->ElevateDegreeForDimension(0);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> original_;
  std::shared_ptr<spl::NURBS<1>> elevated_;
};

TEST_F(ALinearNURBSForDegreeElevationForDimension0, HasELevatedDegreeInDimension0) {  // NOLINT
  ASSERT_THAT(elevated_->GetDegree(0).Get(), original_->GetDegree(0).Get() + 1);
}

TEST_F(ALinearNURBSForDegreeElevationForDimension0, HasNumberOfDifferentKnotsMoreKnots) {  // NOLINT
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(elevated_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots() +
                  original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
}

TEST_F(ALinearNURBSForDegreeElevationForDimension0, HasNumberOfNonZeroKnotSpansMoreControlPoints) {  // NOLINT
  ASSERT_THAT(elevated_->GetTotalNumberOfControlPoints(), original_->GetTotalNumberOfControlPoints() +
      original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
}

TEST_F(ALinearNURBSForDegreeElevationForDimension0, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
}

TEST_F(ALinearNURBSForDegreeElevationForDimension0, WriteRandomSplineToXMLIfAnyPreviousTestFailed) {  // NOLINT
  if (testing::UnitTest::GetInstance()->current_test_case()->failed_test_count() > 0)
    splinelib::test::random_spline_writer::WriteToXML<1>({original_, elevated_}, testing::UnitTest::GetInstance());
}

class A2DBSplineForDegreeElevationForDimension0 : public Test {  // NOLINT
 public:
  A2DBSplineForDegreeElevationForDimension0() {
    std::array<Degree, 2> degree = {Degree{2}, Degree{1}};
    baf::KnotVectors<2> knot_vector_before = {
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.3},
             ParametricCoordinate{0.6},
             ParametricCoordinate{0.6}, ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1}, ParametricCoordinate{1}}))};
    std::vector<spl::ControlPoint> control_points = {
        spl::ControlPoint(std::vector<double>({1.0, 1.0})), spl::ControlPoint(std::vector<double>({1.0, 2.0})),
        spl::ControlPoint(std::vector<double>({2.0, 2.0})), spl::ControlPoint(std::vector<double>({2.0, 3.0})),
        spl::ControlPoint(std::vector<double>({3.0, 3.0})), spl::ControlPoint(std::vector<double>({3.0, 4.0})),

        spl::ControlPoint(std::vector<double>({2.0, 1.0})), spl::ControlPoint(std::vector<double>({1.5, 2.5})),
        spl::ControlPoint(std::vector<double>({2.0, 3.0})), spl::ControlPoint(std::vector<double>({2.5, 3.5})),
        spl::ControlPoint(std::vector<double>({3.5, 4.0})), spl::ControlPoint(std::vector<double>({5.0, 5.0}))
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
  ASSERT_THAT(elevated_->GetDegree(0).Get(), original_->GetDegree(0).Get() + 1);
  ASSERT_THAT(elevated_->GetDegree(1).Get(), original_->GetDegree(1).Get());
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
  ASSERT_THAT(elevated_->GetNumberOfPointsPerDirection()[0],
              original_->GetNumberOfPointsPerDirection()[0] +
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
  ASSERT_THAT(elevated_->GetNumberOfPointsPerDirection()[1], original_->GetNumberOfPointsPerDirection()[1]);
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

TEST_F(A2DBSplineForDegreeElevationForDimension0, WriteRandomSplineToXMLIfAnyPreviousTestFailed) {  // NOLINT
  if (testing::UnitTest::GetInstance()->current_test_case()->failed_test_count() > 0)
    splinelib::test::random_spline_writer::WriteToXML<2>({original_, elevated_}, testing::UnitTest::GetInstance());
}

class Random2DNURBSForDegreeElevationForDimension0 : public Test {  // NOLINT
 public:
  Random2DNURBSForDegreeElevationForDimension0() {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0.0}, ParametricCoordinate{1.0}};
    original_ = splinelib::test::RandomSplineUtils<2>::GenerateRandomNURBS(limits, 3, 3);
    elevated_ = std::make_shared<spl::NURBS<2>>(*original_);
    elevated_->ElevateDegreeForDimension(0);
  }

 protected:
  std::shared_ptr<spl::NURBS<2>> original_;
  std::shared_ptr<spl::NURBS<2>> elevated_;
};

TEST_F(Random2DNURBSForDegreeElevationForDimension0, HasELevatedDegreeInDimension0) {  // NOLINT
  ASSERT_THAT(elevated_->GetDegree(0).Get(), original_->GetDegree(0).Get() + 1);
  ASSERT_THAT(elevated_->GetDegree(1).Get(), original_->GetDegree(1).Get());
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
  ASSERT_THAT(elevated_->GetNumberOfPointsPerDirection()[0],
              original_->GetNumberOfPointsPerDirection()[0] +
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
  ASSERT_THAT(elevated_->GetNumberOfPointsPerDirection()[1], original_->GetNumberOfPointsPerDirection()[1]);
}

TEST_F(Random2DNURBSForDegreeElevationForDimension0, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
}

TEST_F(Random2DNURBSForDegreeElevationForDimension0, WriteRandomSplineToXMLIfAnyPreviousTestFailed) {  // NOLINT
  if (testing::UnitTest::GetInstance()->current_test_case()->failed_test_count() > 0)
    splinelib::test::random_spline_writer::WriteToXML<2>({original_, elevated_}, testing::UnitTest::GetInstance());
}

class Random3DBSplineForDegreeElevationForDimension0 : public Test {  // NOLINT
 public:
  Random3DBSplineForDegreeElevationForDimension0() {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0.0}, ParametricCoordinate{1.0}};
    original_ = splinelib::test::RandomSplineUtils<3>::GenerateRandomBSpline(limits, 4, 3);
    elevated_ = std::make_shared<spl::BSpline<3>>(*original_);
    elevated_->ElevateDegreeForDimension(0);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> original_;
  std::shared_ptr<spl::BSpline<3>> elevated_;
};

TEST_F(Random3DBSplineForDegreeElevationForDimension0, HasELevatedDegreeInDimension0) {  // NOLINT
  ASSERT_THAT(elevated_->GetDegree(0).Get(), original_->GetDegree(0).Get() + 1);
  ASSERT_THAT(elevated_->GetDegree(1).Get(), original_->GetDegree(1).Get());
  ASSERT_THAT(elevated_->GetDegree(2).Get(), original_->GetDegree(2).Get());
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
  ASSERT_THAT(elevated_->GetNumberOfPointsPerDirection()[0],
              original_->GetNumberOfPointsPerDirection()[0] +
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
  ASSERT_THAT(elevated_->GetNumberOfPointsPerDirection()[1], original_->GetNumberOfPointsPerDirection()[1]);
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

TEST_F(Random3DBSplineForDegreeElevationForDimension0, WriteRandomSplineToXMLIfAnyPreviousTestFailed) {  // NOLINT
  if (testing::UnitTest::GetInstance()->current_test_case()->failed_test_count() > 0)
    splinelib::test::random_spline_writer::WriteToXML<3>({original_, elevated_}, testing::UnitTest::GetInstance());
}

class Random3DNURBSForDegreeElevationForDimension1 : public Test {  // NOLINT
 public:
  Random3DNURBSForDegreeElevationForDimension1() {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0.0}, ParametricCoordinate{1.0}};
    original_ = splinelib::test::RandomSplineUtils<3>::GenerateRandomNURBS(limits, 3, 3);
    elevated_ = std::make_shared<spl::NURBS<3>>(*original_);
    elevated_->ElevateDegreeForDimension(1);
  }

 protected:
  std::shared_ptr<spl::NURBS<3>> original_;
  std::shared_ptr<spl::NURBS<3>> elevated_;
};

TEST_F(Random3DNURBSForDegreeElevationForDimension1, HasELevatedDegreeInDimension1) {  // NOLINT
  ASSERT_THAT(elevated_->GetDegree(0).Get(), original_->GetDegree(0).Get());
  ASSERT_THAT(elevated_->GetDegree(1).Get(), original_->GetDegree(1).Get() + 1);
  ASSERT_THAT(elevated_->GetDegree(2).Get(), original_->GetDegree(2).Get());
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
  ASSERT_THAT(elevated_->GetNumberOfPointsPerDirection()[0], original_->GetNumberOfPointsPerDirection()[0]);
  ASSERT_THAT(elevated_->GetNumberOfPointsPerDirection()[1],
              original_->GetNumberOfPointsPerDirection()[1] +
              original_->GetKnotVector(1)->GetNumberOfDifferentKnots() - 1);
}

TEST_F(Random3DNURBSForDegreeElevationForDimension1, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(elevated_->AreGeometricallyEqual(*original_), true);
}

TEST_F(Random3DNURBSForDegreeElevationForDimension1, WriteRandomSplineToXMLIfAnyPreviousTestFailed) {  // NOLINT
  if (testing::UnitTest::GetInstance()->current_test_case()->failed_test_count() > 0)
    splinelib::test::random_spline_writer::WriteToXML<3>({original_, elevated_}, testing::UnitTest::GetInstance());
}
