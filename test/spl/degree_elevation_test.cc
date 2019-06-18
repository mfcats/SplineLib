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

template<int DIM>
void PrintSpline(std::shared_ptr<spl::NURBS<DIM>> spline) {
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
  std::cout << std::endl << "weights:" << std::endl;
  for (int i = 0; i < spline->GetNumberOfControlPoints(); ++i) {
    std::cout << spline->GetWeight({i}) << "  ";
  }
  std::cout << std::endl;
}

using testing::Test;

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

class ALinearNURBS : public Test {  // NOLINT
 public:
  ALinearNURBS() {
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
    std::vector<double> weights = {1.0, 1.5, 0.5, 1.0};
    nurbs_1d_before_ = std::make_shared<spl::NURBS<1>>(knot_vector_before, degree, control_points, weights);
    nurbs_1d_after_ = std::make_shared<spl::NURBS<1>>(knot_vector_after, degree, control_points, weights);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_before_;
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_after_;
};

TEST_F(ALinearNURBS, ElevatesDegreeFrom1To2Correctly) {  // NOLINT
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
  std::vector<double> weights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  auto test = std::make_shared<spl::NURBS<1>>(knot_vector, degree, control_points, weights);

  PrintSpline(nurbs_1d_before_);
  nurbs_1d_after_->ElevateDegree(0);
  PrintSpline(nurbs_1d_after_);
//  ASSERT_THAT(nurbs_1d_after_->AreEqual(*test), true);
//  ASSERT_THAT(nurbs_1d_before_->AreGeometricallyEqual(*test), true);
//  ASSERT_THAT(nurbs_1d_after_->AreGeometricallyEqual(*test), true);
  ASSERT_THAT(nurbs_1d_after_->AreGeometricallyEqual(*nurbs_1d_before_), true);
}

class ALinearBSplineToCompare : public Test {  // NOLINT
 public:
  ALinearBSplineToCompare() {
    std::array<Degree, 1> degree = {Degree{1}};
    KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{1},
                         ParamCoord{1}}))};
    KnotVectors<1> knot_vector_after = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{1},
                         ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.5, 3.0, 1.5})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.5})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0, 1.0}))
    };
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(knot_vector_after, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(ALinearBSplineToCompare, ElevatesDegreeFrom1To2Correctly) {  // NOLINT
  PrintSpline(bspline_1d_before_);
  bspline_1d_after_->ElevateDegree(0);
  PrintSpline(bspline_1d_after_);
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

TEST_F(Random1DBSplineForDegreeElevation, HasELevatedDegree) {  // NOLINT
  ASSERT_THAT(after_elevation->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
}

TEST_F(Random1DBSplineForDegreeElevation, HasMoreKnots) {  // NOLINT
  ASSERT_THAT(after_elevation->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots()
                  + original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
}

TEST_F(Random1DBSplineForDegreeElevation, HasMoreControlPoints) {  // NOLINT
  ASSERT_THAT(after_elevation->GetNumberOfControlPoints(),
              original_->GetNumberOfControlPoints() + original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
}

TEST_F(Random1DBSplineForDegreeElevation, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(after_elevation->AreGeometricallyEqual(*original_), true);
}

class Random1DNURBSForDegreeElevation : public Test {  // NOLINT
 public:
  Random1DNURBSForDegreeElevation() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomNURBSGenerator<1> spline_generator(limits, 5, 3);
    spl::NURBS<1> nurbs(spline_generator);
    original_ = std::make_shared<spl::NURBS<1>>(nurbs);

    spl::NURBS<1> elevation_spline(nurbs);
    elevation_spline.ElevateDegree(0);
    after_elevation_ = std::make_shared<spl::NURBS<1>>(elevation_spline);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> original_;
  std::shared_ptr<spl::NURBS<1>> after_elevation_;
};

TEST_F(Random1DNURBSForDegreeElevation, HasELevatedDegree) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
}

TEST_F(Random1DNURBSForDegreeElevation, HasMoreKnots) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots()
                  + original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
}

TEST_F(Random1DNURBSForDegreeElevation, HasMoreControlPoints) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetNumberOfControlPoints(),
              original_->GetNumberOfControlPoints() + original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
}

TEST_F(Random1DNURBSForDegreeElevation, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
//  io::IRITWriter writer;
//  std::any before = std::make_any<std::shared_ptr<spl::NURBS<1>>>(original_);
//  std::any after = std::make_any<std::shared_ptr<spl::NURBS<1>>>(after_elevation_);
//  writer.WriteFile({before, after}, "degree_elevation.itd");
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_, 0.01), true);
}

class A2DBSplineForDegreeElevation : public Test {  // NOLINT
 public:
  A2DBSplineForDegreeElevation() {
    std::array<Degree, 2> degree = {Degree{2}, Degree{1}};
    KnotVectors<2> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{0.6},
                         ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
                                         std::make_shared<baf::KnotVector>(
                                             baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                              ParamCoord{1}}))};
    KnotVectors<2> knot_vector_after = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{0.6},
                         ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
                                        std::make_shared<baf::KnotVector>(
                                            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                             ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0})),
        baf::ControlPoint(std::vector<double>({3.0, 3.0})),
        baf::ControlPoint(std::vector<double>({3.0, 4.0})),

        baf::ControlPoint(std::vector<double>({2.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.5, 2.5})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0})),
        baf::ControlPoint(std::vector<double>({2.5, 3.5})),
        baf::ControlPoint(std::vector<double>({3.5, 4.0})),
        baf::ControlPoint(std::vector<double>({5.0, 5.0}))
    };
    original_ = std::make_shared<spl::BSpline<2>>(knot_vector_before, degree, control_points);
    after_elevation_ = std::make_shared<spl::BSpline<2>>(knot_vector_after, degree, control_points);
    after_elevation_->ElevateDegree(0);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> original_;
  std::shared_ptr<spl::BSpline<2>> after_elevation_;
};

TEST_F(A2DBSplineForDegreeElevation, HasELevatedDegree) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get());
}

TEST_F(A2DBSplineForDegreeElevation, HasMoreKnots) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots()
                  + original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(), original_->GetKnotVector(1)->GetNumberOfKnots());
}

TEST_F(A2DBSplineForDegreeElevation, HasMoreControlPoints) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0],
              original_->GetPointsPerDirection()[0] + original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
}

TEST_F(A2DBSplineForDegreeElevation, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

TEST_F(A2DBSplineForDegreeElevation, HasELevatedDegreeInDirection0and1) {  // NOLINT
  after_elevation_->ElevateDegree(1);
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get() + 1);

  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots()
                  + original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots()
                  + original_->GetKnotVector(1)->GetNumberOfDifferentKnots());

  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0],
              original_->GetPointsPerDirection()[0] + original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1],
              original_->GetPointsPerDirection()[1] + original_->GetKnotVector(1)->GetNumberOfDifferentKnots() - 1);

  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

TEST_F(A2DBSplineForDegreeElevation, HasELevatedDegreeTwiceInDirection0) {  // NOLINT
  after_elevation_->ElevateDegree(0);
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get() + 2);
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get());

  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots()
                  + 2 * original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(), original_->GetKnotVector(1)->GetNumberOfKnots());

  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0],
              original_->GetPointsPerDirection()[0]
                  + 2 * (original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1));
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);

  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

TEST_F(A2DBSplineForDegreeElevation, HasELevatedDegreeTwiceInDirection0andOnceInDirection1) {  // NOLINT
  after_elevation_->ElevateDegree(0);
  after_elevation_->ElevateDegree(1);
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get() + 2);
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get() + 1);

  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots()
                  + 2 * original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots()
                  + 1 * original_->GetKnotVector(1)->GetNumberOfDifferentKnots());

  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0],
              original_->GetPointsPerDirection()[0]
                  + 2 * (original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1));
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1],
              original_->GetPointsPerDirection()[1]
                  + 1 * (original_->GetKnotVector(1)->GetNumberOfDifferentKnots() - 1));

  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_, 1e-10), true);
}

TEST_F(A2DBSplineForDegreeElevation, HasELevatedDegreeTwiceInDirection0and1) {  // NOLINT
  after_elevation_->ElevateDegree(0);
  after_elevation_->ElevateDegree(1);
  after_elevation_->ElevateDegree(1);
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get() + 2);
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get() + 2);

  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots()
                  + 2 * original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots()
                  + 2 * original_->GetKnotVector(1)->GetNumberOfDifferentKnots());

  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0],
              original_->GetPointsPerDirection()[0]
                  + 2 * (original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1));
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1],
              original_->GetPointsPerDirection()[1]
                  + 2 * (original_->GetKnotVector(1)->GetNumberOfDifferentKnots() - 1));

  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_, 1e-12), true);
}

class A2DBSplineForDegreeElevationInDirection1 : public Test {  // NOLINT
 public:
  A2DBSplineForDegreeElevationInDirection1() {
    std::array<Degree, 2> degree = {Degree{1}, Degree{2}};
    KnotVectors<2> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                         ParamCoord{1}})),
                                         std::make_shared<baf::KnotVector>(
                                             baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                              ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{0.6},
                                                              ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    KnotVectors<2> knot_vector_after = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                         ParamCoord{1}})),
                                        std::make_shared<baf::KnotVector>(
                                            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                             ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{0.6},
                                                             ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0})),

        baf::ControlPoint(std::vector<double>({1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.5, 2.5})),

        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0})),

        baf::ControlPoint(std::vector<double>({2.0, 3.0})),
        baf::ControlPoint(std::vector<double>({2.5, 3.5})),

        baf::ControlPoint(std::vector<double>({3.0, 3.0})),
        baf::ControlPoint(std::vector<double>({3.5, 4.0})),

        baf::ControlPoint(std::vector<double>({3.0, 4.0})),
        baf::ControlPoint(std::vector<double>({5.0, 5.0}))
    };
    original_ = std::make_shared<spl::BSpline<2>>(knot_vector_before, degree, control_points);
    after_elevation_ = std::make_shared<spl::BSpline<2>>(knot_vector_after, degree, control_points);
    after_elevation_->ElevateDegree(1);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> original_;
  std::shared_ptr<spl::BSpline<2>> after_elevation_;
};

TEST_F(A2DBSplineForDegreeElevationInDirection1, HasELevatedDegree) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get());
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get() + 1);
}

TEST_F(A2DBSplineForDegreeElevationInDirection1, HasMoreKnots) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(1)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(), original_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots()
                  + original_->GetKnotVector(1)->GetNumberOfDifferentKnots());
}

TEST_F(A2DBSplineForDegreeElevationInDirection1, HasMoreControlPoints) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0]);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1],
              original_->GetPointsPerDirection()[1] + original_->GetKnotVector(1)->GetNumberOfDifferentKnots() - 1);
}

TEST_F(A2DBSplineForDegreeElevationInDirection1, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

class Random2DBSplineForDegreeElevationInDirection0 : public Test {  // NOLINT
 public:
  Random2DBSplineForDegreeElevationInDirection0() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomBSplineGenerator<2> spline_generator(limits, 2, 3);
    spl::BSpline<2> b_spline(spline_generator);
    original_ = std::make_shared<spl::BSpline<2>>(b_spline);

    spl::BSpline<2> elevation_spline(b_spline);
    elevation_spline.ElevateDegree(0);
    after_elevation_ = std::make_shared<spl::BSpline<2>>(elevation_spline);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> original_;
  std::shared_ptr<spl::BSpline<2>> after_elevation_;
};

TEST_F(Random2DBSplineForDegreeElevationInDirection0, HasELevatedDegree) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get());
}

TEST_F(Random2DBSplineForDegreeElevationInDirection0, HasMoreKnots) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots()
                  + original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(), original_->GetKnotVector(1)->GetNumberOfKnots());
}

TEST_F(Random2DBSplineForDegreeElevationInDirection0, HasMoreControlPoints) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0],
              original_->GetPointsPerDirection()[0] + original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
}

TEST_F(Random2DBSplineForDegreeElevationInDirection0, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

TEST_F(Random2DBSplineForDegreeElevationInDirection0, ElevatesDegreeTwiceInDirection0) {  // NOLINT
  after_elevation_->ElevateDegree(0);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

TEST_F(Random2DBSplineForDegreeElevationInDirection0, ElevatesDegreeInDirection0and1) {  // NOLINT
  after_elevation_->ElevateDegree(1);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

TEST_F(Random2DBSplineForDegreeElevationInDirection0, ElevatesDegreeInDirection0andTwiceInDirection1) {  // NOLINT
  after_elevation_->ElevateDegree(1);
  after_elevation_->ElevateDegree(1);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

TEST_F(Random2DBSplineForDegreeElevationInDirection0, ElevatesDegreeTwiceInDirection0and1) {  // NOLINT
  after_elevation_->ElevateDegree(0);
  after_elevation_->ElevateDegree(1);
  after_elevation_->ElevateDegree(1);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

class Random2DBSplineForDegreeElevationInDirection1 : public Test {  // NOLINT
 public:
  Random2DBSplineForDegreeElevationInDirection1() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomBSplineGenerator<2> spline_generator(limits, 2, 3);
    spl::BSpline<2> b_spline(spline_generator);
    original_ = std::make_shared<spl::BSpline<2>>(b_spline);

    spl::BSpline<2> elevation_spline(b_spline);
    elevation_spline.ElevateDegree(1);
    after_elevation_ = std::make_shared<spl::BSpline<2>>(elevation_spline);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> original_;
  std::shared_ptr<spl::BSpline<2>> after_elevation_;
};

TEST_F(Random2DBSplineForDegreeElevationInDirection1, HasELevatedDegree) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get());
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get() + 1);
}

TEST_F(Random2DBSplineForDegreeElevationInDirection1, HasMoreKnots) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(), original_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots()
                  + original_->GetKnotVector(1)->GetNumberOfDifferentKnots());
}

TEST_F(Random2DBSplineForDegreeElevationInDirection1, HasMoreControlPoints) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0]);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1],
              original_->GetPointsPerDirection()[1] + original_->GetKnotVector(1)->GetNumberOfDifferentKnots() - 1);
}

TEST_F(Random2DBSplineForDegreeElevationInDirection1, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

class A3DBSplineForDegreeElevation : public Test {  // NOLINT
 public:
  A3DBSplineForDegreeElevation() {
    std::array<Degree, 3> degree = {Degree{2}, Degree{1}, Degree{1}};
    KnotVectors<3> knot_vector_before = {
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{0.6},
             ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}}))};
    KnotVectors<3> knot_vector_after = {
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6},
                             ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0, 0.0})),
        baf::ControlPoint(std::vector<double>({3.0, 3.0, 0.0})),
        baf::ControlPoint(std::vector<double>({3.0, 4.0, 0.0})),

        baf::ControlPoint(std::vector<double>({2.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.5, 2.5, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.5, 3.5, 0.0})),
        baf::ControlPoint(std::vector<double>({3.5, 4.0, 0.0})),
        baf::ControlPoint(std::vector<double>({5.0, 5.0, 0.0})),

        baf::ControlPoint(std::vector<double>({1.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.0, 1.5})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0, 2.5})),
        baf::ControlPoint(std::vector<double>({3.0, 3.0, 1.5})),
        baf::ControlPoint(std::vector<double>({3.0, 4.0, 1.0})),

        baf::ControlPoint(std::vector<double>({2.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.5, 2.5, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0, 2.5})),
        baf::ControlPoint(std::vector<double>({2.5, 3.5, 3.0})),
        baf::ControlPoint(std::vector<double>({3.5, 4.0, 2.0})),
        baf::ControlPoint(std::vector<double>({5.0, 5.0, 1.0}))
    };
    original_ = std::make_shared<spl::BSpline<3>>(knot_vector_before, degree, control_points);
    after_elevation_ = std::make_shared<spl::BSpline<3>>(knot_vector_after, degree, control_points);
    after_elevation_->ElevateDegree(0);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> original_;
  std::shared_ptr<spl::BSpline<3>> after_elevation_;
};

TEST_F(A3DBSplineForDegreeElevation, HasELevatedDegree) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get());
  ASSERT_THAT(after_elevation_->GetDegree(2).get(), original_->GetDegree(2).get());
}

TEST_F(A3DBSplineForDegreeElevation, HasMoreKnots) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots()
                  + original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(), original_->GetKnotVector(1)->GetNumberOfKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(2)->GetNumberOfKnots(), original_->GetKnotVector(2)->GetNumberOfKnots());
}

TEST_F(A3DBSplineForDegreeElevation, HasMoreControlPoints) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0],
              original_->GetPointsPerDirection()[0] + original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[2], original_->GetPointsPerDirection()[2]);
}

TEST_F(A3DBSplineForDegreeElevation, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_, 1e-12), true);
}

TEST_F(A3DBSplineForDegreeElevation, ElevatesDegreeInDirection0and1and2) {  // NOLINT
  after_elevation_->ElevateDegree(1);
  after_elevation_->ElevateDegree(2);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_, 1e-12), true);
}

TEST_F(A3DBSplineForDegreeElevation, ElevatesDegreeInDirection0and1) {  // NOLINT
  after_elevation_->ElevateDegree(1);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_, 1e-12), true);
}

TEST_F(A3DBSplineForDegreeElevation, ElevatesDegreeInDirection0and2) {  // NOLINT
  after_elevation_->ElevateDegree(2);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_, 1e-12), true);
}

class A3DBSplineForDegreeElevationInDirection1 : public Test {  // NOLINT
 public:
  A3DBSplineForDegreeElevationInDirection1() {
    std::array<Degree, 3> degree = {Degree{1}, Degree{2}, Degree{1}};
    KnotVectors<3> knot_vector_before = {
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{0.6},
             ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}}))};
    KnotVectors<3> knot_vector_after = {
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6},
                             ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0, 0.0})),

        baf::ControlPoint(std::vector<double>({1.0, 2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.5, 2.5, 0.0})),

        baf::ControlPoint(std::vector<double>({2.0, 2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0, 0.0})),

        baf::ControlPoint(std::vector<double>({2.0, 3.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.5, 3.5, 0.0})),

        baf::ControlPoint(std::vector<double>({3.0, 3.0, 0.0})),
        baf::ControlPoint(std::vector<double>({3.5, 4.0, 0.0})),

        baf::ControlPoint(std::vector<double>({3.0, 4.0, 0.0})),
        baf::ControlPoint(std::vector<double>({5.0, 5.0, 0.0})),

        baf::ControlPoint(std::vector<double>({1.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0, 1.0})),

        baf::ControlPoint(std::vector<double>({1.0, 2.0, 1.5})),
        baf::ControlPoint(std::vector<double>({1.5, 2.5, 2.0})),

        baf::ControlPoint(std::vector<double>({2.0, 2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0, 2.5})),

        baf::ControlPoint(std::vector<double>({2.0, 3.0, 2.5})),
        baf::ControlPoint(std::vector<double>({2.5, 3.5, 3.0})),

        baf::ControlPoint(std::vector<double>({3.0, 3.0, 1.5})),
        baf::ControlPoint(std::vector<double>({3.5, 4.0, 2.0})),

        baf::ControlPoint(std::vector<double>({3.0, 4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({5.0, 5.0, 1.0}))
    };
    original_ = std::make_shared<spl::BSpline<3>>(knot_vector_before, degree, control_points);
    after_elevation_ = std::make_shared<spl::BSpline<3>>(knot_vector_after, degree, control_points);
    after_elevation_->ElevateDegree(1);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> original_;
  std::shared_ptr<spl::BSpline<3>> after_elevation_;
};

TEST_F(A3DBSplineForDegreeElevationInDirection1, HasELevatedDegree) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get());
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get() + 1);
  ASSERT_THAT(after_elevation_->GetDegree(2).get(), original_->GetDegree(2).get());
}

TEST_F(A3DBSplineForDegreeElevationInDirection1, HasMoreKnots) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(), original_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots()
                  + original_->GetKnotVector(1)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(2)->GetNumberOfKnots(), original_->GetKnotVector(2)->GetNumberOfKnots());
}

TEST_F(A3DBSplineForDegreeElevationInDirection1, HasMoreControlPoints) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0]);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1],
              original_->GetPointsPerDirection()[1] + original_->GetKnotVector(1)->GetNumberOfDifferentKnots() - 1);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[2], original_->GetPointsPerDirection()[2]);
}

TEST_F(A3DBSplineForDegreeElevationInDirection1, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

TEST_F(A3DBSplineForDegreeElevationInDirection1, ElevatesDegreeInDirection0and1) {  // NOLINT
  after_elevation_->ElevateDegree(0);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_, 1e-10), true);
}

TEST_F(A3DBSplineForDegreeElevationInDirection1, ElevatesDegreeInDirection1and2) {  // NOLINT
  after_elevation_->ElevateDegree(2);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_, 1e-10), true);
}

TEST_F(A3DBSplineForDegreeElevationInDirection1, ElevatesDegreeInDirection0and1and2) {  // NOLINT
  after_elevation_->ElevateDegree(0);
  after_elevation_->ElevateDegree(2);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_, 1e-10), true);
}

class A3DBSplineForDegreeElevationInDirection2 : public Test {  // NOLINT
 public:
  A3DBSplineForDegreeElevationInDirection2() {
    std::array<Degree, 3> degree = {Degree{1}, Degree{1}, Degree{2}};
    KnotVectors<3> knot_vector_before = {
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector(
            {ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{0.6},
             ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    KnotVectors<3> knot_vector_after = {
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6},
                             ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0, 0.0})),

        baf::ControlPoint(std::vector<double>({1.0, 1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0, 1.0})),

        baf::ControlPoint(std::vector<double>({1.0, 2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.5, 2.5, 0.0})),

        baf::ControlPoint(std::vector<double>({1.0, 2.0, 1.5})),
        baf::ControlPoint(std::vector<double>({1.5, 2.5, 2.0})),

        baf::ControlPoint(std::vector<double>({2.0, 2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0, 0.0})),

        baf::ControlPoint(std::vector<double>({2.0, 2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 3.0, 2.5})),

        baf::ControlPoint(std::vector<double>({2.0, 3.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.5, 3.5, 0.0})),

        baf::ControlPoint(std::vector<double>({2.0, 3.0, 2.5})),
        baf::ControlPoint(std::vector<double>({2.5, 3.5, 3.0})),

        baf::ControlPoint(std::vector<double>({3.0, 3.0, 0.0})),
        baf::ControlPoint(std::vector<double>({3.5, 4.0, 0.0})),

        baf::ControlPoint(std::vector<double>({3.0, 3.0, 1.5})),
        baf::ControlPoint(std::vector<double>({3.5, 4.0, 2.0})),

        baf::ControlPoint(std::vector<double>({3.0, 4.0, 0.0})),
        baf::ControlPoint(std::vector<double>({5.0, 5.0, 0.0})),

        baf::ControlPoint(std::vector<double>({3.0, 4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({5.0, 5.0, 1.0}))
    };
    original_ = std::make_shared<spl::BSpline<3>>(knot_vector_before, degree, control_points);
    after_elevation_ = std::make_shared<spl::BSpline<3>>(knot_vector_after, degree, control_points);
    after_elevation_->ElevateDegree(2);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> original_;
  std::shared_ptr<spl::BSpline<3>> after_elevation_;
};

TEST_F(A3DBSplineForDegreeElevationInDirection2, HasELevatedDegree) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get());
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get());
  ASSERT_THAT(after_elevation_->GetDegree(2).get(), original_->GetDegree(2).get() + 1);
}

TEST_F(A3DBSplineForDegreeElevationInDirection2, HasMoreKnots) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(), original_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(), original_->GetKnotVector(1)->GetNumberOfKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(2)->GetNumberOfKnots(),
              original_->GetKnotVector(2)->GetNumberOfKnots()
                  + original_->GetKnotVector(2)->GetNumberOfDifferentKnots());
}

TEST_F(A3DBSplineForDegreeElevationInDirection2, HasMoreControlPoints) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0]);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[2],
              original_->GetPointsPerDirection()[2] + original_->GetKnotVector(2)->GetNumberOfDifferentKnots() - 1);
}

TEST_F(A3DBSplineForDegreeElevationInDirection2, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

TEST_F(A3DBSplineForDegreeElevationInDirection2, ElevatesDegreeInDirection0and2) {  // NOLINT
  after_elevation_->ElevateDegree(0);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_, 1e-12), true);
}

TEST_F(A3DBSplineForDegreeElevationInDirection2, ElevatesDegreeInDirection1and2) {  // NOLINT
  after_elevation_->ElevateDegree(1);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_, 1e-12), true);
}

TEST_F(A3DBSplineForDegreeElevationInDirection2, ElevatesDegreeInDirection0and1and2) {  // NOLINT
  after_elevation_->ElevateDegree(0);
  after_elevation_->ElevateDegree(1);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_, 1e-12), true);
}

class Random3DBSplineForDegreeElevationInDirection0 : public Test {  // NOLINT
 public:
  Random3DBSplineForDegreeElevationInDirection0() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomBSplineGenerator<3> spline_generator(limits, 4, 3);
    spl::BSpline<3> b_spline(spline_generator);
    original_ = std::make_shared<spl::BSpline<3>>(b_spline);

    spl::BSpline<3> elevation_spline(b_spline);
    elevation_spline.ElevateDegree(0);
    after_elevation_ = std::make_shared<spl::BSpline<3>>(elevation_spline);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> original_;
  std::shared_ptr<spl::BSpline<3>> after_elevation_;
};

TEST_F(Random3DBSplineForDegreeElevationInDirection0, HasELevatedDegree) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get());
  ASSERT_THAT(after_elevation_->GetDegree(2).get(), original_->GetDegree(2).get());
}

TEST_F(Random3DBSplineForDegreeElevationInDirection0, HasMoreKnots) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots()
                  + original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(), original_->GetKnotVector(1)->GetNumberOfKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(2)->GetNumberOfKnots(), original_->GetKnotVector(2)->GetNumberOfKnots());
}

TEST_F(Random3DBSplineForDegreeElevationInDirection0, HasMoreControlPoints) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0],
              original_->GetPointsPerDirection()[0] + original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
}

TEST_F(Random3DBSplineForDegreeElevationInDirection0, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

TEST_F(Random3DBSplineForDegreeElevationInDirection0, ElevatesDegreeInAllDirections) {  // NOLINT
  after_elevation_->ElevateDegree(1);
  after_elevation_->ElevateDegree(2);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

class Random3DNURBSForDegreeElevationInDirection0 : public Test {  // NOLINT
 public:
  Random3DNURBSForDegreeElevationInDirection0() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomNURBSGenerator<3> spline_generator(limits, 4, 3);
    spl::NURBS<3> b_spline(spline_generator);
    original_ = std::make_shared<spl::NURBS<3>>(b_spline);

    spl::NURBS<3> elevation_spline(b_spline);
    elevation_spline.ElevateDegree(0);
    after_elevation_ = std::make_shared<spl::NURBS<3>>(elevation_spline);
  }

 protected:
  std::shared_ptr<spl::NURBS<3>> original_;
  std::shared_ptr<spl::NURBS<3>> after_elevation_;
};

TEST_F(Random3DNURBSForDegreeElevationInDirection0, HasELevatedDegree) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get() + 1);
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get());
  ASSERT_THAT(after_elevation_->GetDegree(2).get(), original_->GetDegree(2).get());
}

TEST_F(Random3DNURBSForDegreeElevationInDirection0, HasMoreKnots) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots()
                  + original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(), original_->GetKnotVector(1)->GetNumberOfKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(2)->GetNumberOfKnots(), original_->GetKnotVector(2)->GetNumberOfKnots());
}

TEST_F(Random3DNURBSForDegreeElevationInDirection0, HasMoreControlPoints) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0],
              original_->GetPointsPerDirection()[0] + original_->GetKnotVector(0)->GetNumberOfDifferentKnots() - 1);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
}

TEST_F(Random3DNURBSForDegreeElevationInDirection0, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

TEST_F(Random3DNURBSForDegreeElevationInDirection0, ElevatesDegreeInAllDirections) {  // NOLINT
  after_elevation_->ElevateDegree(1);
  after_elevation_->ElevateDegree(2);
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

class Random4DBSplineForDegreeElevationInDirection3 : public Test {  // NOLINT
 public:
  Random4DBSplineForDegreeElevationInDirection3() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomBSplineGenerator<4> spline_generator(limits, 2, 3);
    spl::BSpline<4> b_spline(spline_generator);
    original_ = std::make_shared<spl::BSpline<4>>(b_spline);

    spl::BSpline<4> elevation_spline(b_spline);
    elevation_spline.ElevateDegree(3);
    after_elevation_ = std::make_shared<spl::BSpline<4>>(elevation_spline);
  }

 protected:
  std::shared_ptr<spl::BSpline<4>> original_;
  std::shared_ptr<spl::BSpline<4>> after_elevation_;
};

TEST_F(Random4DBSplineForDegreeElevationInDirection3, HasELevatedDegree) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetDegree(0).get(), original_->GetDegree(0).get());
  ASSERT_THAT(after_elevation_->GetDegree(1).get(), original_->GetDegree(1).get());
  ASSERT_THAT(after_elevation_->GetDegree(2).get(), original_->GetDegree(2).get());
  ASSERT_THAT(after_elevation_->GetDegree(3).get(), original_->GetDegree(3).get() + 1);
}

TEST_F(Random4DBSplineForDegreeElevationInDirection3, HasMoreKnots) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfDifferentKnots(),
              original_->GetKnotVector(0)->GetNumberOfDifferentKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(0)->GetNumberOfKnots(), original_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(1)->GetNumberOfKnots(), original_->GetKnotVector(1)->GetNumberOfKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(2)->GetNumberOfKnots(), original_->GetKnotVector(2)->GetNumberOfKnots());
  ASSERT_THAT(after_elevation_->GetKnotVector(3)->GetNumberOfKnots(),
              original_->GetKnotVector(3)->GetNumberOfKnots()
                  + original_->GetKnotVector(3)->GetNumberOfDifferentKnots());
}

TEST_F(Random4DBSplineForDegreeElevationInDirection3, HasMoreControlPoints) {  // NOLINT
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0]);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1]);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[2], original_->GetPointsPerDirection()[2]);
  ASSERT_THAT(after_elevation_->GetPointsPerDirection()[3],
              original_->GetPointsPerDirection()[3] + original_->GetKnotVector(3)->GetNumberOfDifferentKnots() - 1);
}

TEST_F(Random4DBSplineForDegreeElevationInDirection3, DoesNotChangeGeometricallyAfterDegreeElevation) {  // NOLINT
  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
}

//TEST_F(Random4DBSplineForDegreeElevationInDirection0, ElevatesDegreeInAllDirections) {  // NOLINT
//  after_elevation_->ElevateDegree(1);
//  after_elevation_->ElevateDegree(2);
//  ASSERT_THAT(after_elevation_->AreGeometricallyEqual(*original_), true);
//}
