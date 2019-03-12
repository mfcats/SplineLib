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
using testing::DoubleNear;

class Random2DBSplineForKnotInsertion : public Test {  // NOLINT
 public:
  Random2DBSplineForKnotInsertion() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomBSplineGenerator<2> spline_generator(limits, 10, 3);
    spl::BSpline<2> b_spline(spline_generator);
    bspline_2d_before_ = std::make_shared<spl::BSpline<2>>(b_spline);
    spl::BSpline<2> b_spline_after(b_spline);
    bspline_2d_after_ = std::make_shared<spl::BSpline<2>>(b_spline_after);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> bspline_2d_before_;
  std::shared_ptr<spl::BSpline<2>> bspline_2d_after_;
};

TEST_F(Random2DBSplineForKnotInsertion, InsertsKnot0_4InFirstDirectionCorrectly) {  // NOLINT
  bspline_2d_after_->InsertKnot(ParamCoord(0.4), 0);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(1)->GetNumberOfKnots());
  ASSERT_THAT(bspline_2d_after_->GetPointsPerDirection()[0], bspline_2d_before_->GetPointsPerDirection()[0] + 1);
  double steps = 50;
  for (int i = 0; i <= steps; ++i) {
    ParamCoord coord2 = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    std::array<ParamCoord, 2> param_coord{ParamCoord(i / steps), coord2};
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(bspline_2d_after_->Evaluate(param_coord, {j})[0],
                  DoubleNear(bspline_2d_before_->Evaluate(param_coord, {j})[0], 0.000001));
    }
  }
}

TEST_F(Random2DBSplineForKnotInsertion, InsertsKnot0_7InSecondDirectionCorrectly) {  // NOLINT
  bspline_2d_after_->InsertKnot(ParamCoord(0.7), 1);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(1)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_2d_after_->GetPointsPerDirection()[1], bspline_2d_before_->GetPointsPerDirection()[1] + 1);
  double steps = 50;
  for (int i = 0; i <= steps; ++i) {
    ParamCoord coord1 = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    std::array<ParamCoord, 2> param_coord{coord1, ParamCoord(i / steps)};
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(bspline_2d_after_->Evaluate(param_coord, {j})[0],
                  DoubleNear(bspline_2d_before_->Evaluate(param_coord, {j})[0], 0.000001));
    }
  }
}

TEST_F(Random2DBSplineForKnotInsertion, InsertsKnot0_4InFirstAndKnot0_7InSecondDirectionCorrectly) {  // NOLINT
  bspline_2d_after_->InsertKnot(ParamCoord(0.4), 0);
  bspline_2d_after_->InsertKnot(ParamCoord(0.7), 1);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_2d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              bspline_2d_before_->GetKnotVector(1)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_2d_after_->GetPointsPerDirection()[0], bspline_2d_before_->GetPointsPerDirection()[0] + 1);
  ASSERT_THAT(bspline_2d_after_->GetPointsPerDirection()[1], bspline_2d_before_->GetPointsPerDirection()[1] + 1);
  double steps = 10;
  for (int i = 0; i <= steps; ++i) {
    for (int j = 0; j <= steps; ++j) {
      std::array<ParamCoord, 2> param_coord{ParamCoord(i / steps), ParamCoord(j / steps)};
      for (int k = 0; k < 2; ++k) {
        ASSERT_THAT(bspline_2d_after_->Evaluate(param_coord, {k})[0],
                    DoubleNear(bspline_2d_before_->Evaluate(param_coord, {k})[0], 0.000001));
      }
    }
  }
}

class Random3DBSplineForKnotInsertion : public Test {  // NOLINT
 public:
  Random3DBSplineForKnotInsertion() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomBSplineGenerator<3> spline_generator(limits, 10, 4);
    spl::BSpline<3> b_spline(spline_generator);
    bspline_3d_before_ = std::make_shared<spl::BSpline<3>>(b_spline);
    spl::BSpline<3> b_spline_after(b_spline);
    bspline_3d_after_ = std::make_shared<spl::BSpline<3>>(b_spline_after);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> bspline_3d_before_;
  std::shared_ptr<spl::BSpline<3>> bspline_3d_after_;
};

TEST_F(Random3DBSplineForKnotInsertion, InsertsKnot0_4InFirst_Knot0_99InSecond_Knot0_01InThirdDirection) {  // NOLINT
  bspline_3d_after_->InsertKnot(ParamCoord(0.4), 0);
  bspline_3d_after_->InsertKnot(ParamCoord(0.99), 1);
  bspline_3d_after_->InsertKnot(ParamCoord(0.01), 2);
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              bspline_3d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              bspline_3d_before_->GetKnotVector(1)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_3d_after_->GetKnotVector(2)->GetNumberOfKnots(),
              bspline_3d_before_->GetKnotVector(2)->GetNumberOfKnots() + 1);
  ASSERT_THAT(bspline_3d_after_->GetPointsPerDirection()[0], bspline_3d_before_->GetPointsPerDirection()[0] + 1);
  ASSERT_THAT(bspline_3d_after_->GetPointsPerDirection()[1], bspline_3d_before_->GetPointsPerDirection()[1] + 1);
  ASSERT_THAT(bspline_3d_after_->GetPointsPerDirection()[2], bspline_3d_before_->GetPointsPerDirection()[2] + 1);
  double s = 5;
  for (int i = 0; i <= s; ++i) {
    for (int j = 0; j <= s; ++j) {
      for (int k = 0; k <= s; ++k) {
        std::array<ParamCoord, 3> param_coord{ParamCoord(i / s), ParamCoord(j / s), ParamCoord(k / s)};
        for (int l = 0; l < 3; ++l) {
          ASSERT_THAT(bspline_3d_after_->Evaluate(param_coord, {l})[0],
                      DoubleNear(bspline_3d_before_->Evaluate(param_coord, {l})[0], 0.000001));
        }
      }
    }
  }
}

class Random3DNURBSForKnotInsertion : public Test {  // NOLINT
 public:
  Random3DNURBSForKnotInsertion() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomNURBSGenerator<3> spline_generator(limits, 10, 4);
    spl::NURBS<3> nurbs(spline_generator);
    nurbs_3d_before_ = std::make_shared<spl::NURBS<3>>(nurbs);
    spl::NURBS<3> nurbs_after(nurbs);
    nurbs_3d_after_ = std::make_shared<spl::NURBS<3>>(nurbs_after);
  }

 protected:
  std::shared_ptr<spl::NURBS<3>> nurbs_3d_before_;
  std::shared_ptr<spl::NURBS<3>> nurbs_3d_after_;
};

TEST_F(Random3DNURBSForKnotInsertion, InsertsKnot0_4InFirst_Knot0_99InSecond_Knot0_01InThirdDirection) {  // NOLINT
  nurbs_3d_after_->InsertKnot(ParamCoord(0.4), 0);
  nurbs_3d_after_->InsertKnot(ParamCoord(0.99), 1);
  nurbs_3d_after_->InsertKnot(ParamCoord(0.01), 2);
  ASSERT_THAT(nurbs_3d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_3d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(nurbs_3d_after_->GetKnotVector(1)->GetNumberOfKnots(),
              nurbs_3d_before_->GetKnotVector(1)->GetNumberOfKnots() + 1);
  ASSERT_THAT(nurbs_3d_after_->GetKnotVector(2)->GetNumberOfKnots(),
              nurbs_3d_before_->GetKnotVector(2)->GetNumberOfKnots() + 1);
  ASSERT_THAT(nurbs_3d_after_->GetWeights().size(), nurbs_3d_after_->GetNumberOfControlPoints());
  ASSERT_THAT(nurbs_3d_after_->GetPointsPerDirection()[0], nurbs_3d_before_->GetPointsPerDirection()[0] + 1);
  ASSERT_THAT(nurbs_3d_after_->GetPointsPerDirection()[1], nurbs_3d_before_->GetPointsPerDirection()[1] + 1);
  ASSERT_THAT(nurbs_3d_after_->GetPointsPerDirection()[2], nurbs_3d_before_->GetPointsPerDirection()[2] + 1);
  std::array<ParamCoord, 3> param_coord{};
  for (int i = 0; i < 3; ++i) {
    param_coord[i] = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
  }
  for (int l = 0; l < 3; ++l) {
    ASSERT_THAT(nurbs_3d_after_->Evaluate(param_coord, {l})[0],
                DoubleNear(nurbs_3d_before_->Evaluate(param_coord, {l})[0], 0.000001));
  }
}

class Random1DNURBSForKnotRemoval : public Test {  // NOLINT
 public:
  Random1DNURBSForKnotRemoval() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomNURBSGenerator<1> spline_generator(limits, 10, 3);
    spl::NURBS<1> nurbs(spline_generator);
    nurbs_1d_before_ = std::make_shared<spl::NURBS<1>>(nurbs);
    spl::NURBS<1> nurbs_after(nurbs);
    nurbs_1d_after_ = std::make_shared<spl::NURBS<1>>(nurbs_after);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_before_;
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_after_;
};

TEST_F(Random1DNURBSForKnotRemoval, InsertsAndRemovesKnot0_5) {  // NOLINT
  nurbs_1d_after_->InsertKnot(ParamCoord(0.5), 0);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_1d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  nurbs_1d_after_->RemoveKnot(ParamCoord(0.5), 0, 1e-15);
  ASSERT_THAT(nurbs_1d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_1d_before_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(nurbs_1d_after_->GetNumberOfControlPoints(), nurbs_1d_before_->GetNumberOfControlPoints());

  for (int point = 0; point < nurbs_1d_before_->GetNumberOfControlPoints(); ++point) {
    for (int dimension = 0; dimension < nurbs_1d_before_->GetDimension(); ++dimension) {
      ASSERT_THAT(nurbs_1d_after_->GetControlPoint({point}, dimension),
                  DoubleNear(nurbs_1d_before_->GetControlPoint({point}, dimension), 1e-12));
    }
    ASSERT_THAT(nurbs_1d_after_->GetWeight({point}), DoubleNear(nurbs_1d_before_->GetWeight({point}), 1e-12));
  }
  std::array<ParamCoord, 1> param_coord{};
  for (int evaluated_point = 0; evaluated_point < 100; ++evaluated_point) {
    param_coord[0] = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    for (int dimension = 0; dimension < nurbs_1d_before_->GetDimension(); ++dimension) {
      ASSERT_THAT(nurbs_1d_after_->Evaluate(param_coord, {dimension})[0],
                  DoubleNear(nurbs_1d_before_->Evaluate(param_coord, {dimension})[0], 1e-12));
    }
  }
}

class Random2DNURBSForKnotRemoval : public Test {  // NOLINT
 public:
  Random2DNURBSForKnotRemoval() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomNURBSGenerator<2> spline_generator(limits, 8, 3);
    spl::NURBS<2> nurbs(spline_generator);
    nurbs_2d_before_ = std::make_shared<spl::NURBS<2>>(nurbs);
    spl::NURBS<2> nurbs_after(nurbs);
    nurbs_2d_after_ = std::make_shared<spl::NURBS<2>>(nurbs_after);
  }

 protected:
  std::shared_ptr<spl::NURBS<2>> nurbs_2d_before_;
  std::shared_ptr<spl::NURBS<2>> nurbs_2d_after_;
};

TEST_F(Random2DNURBSForKnotRemoval, InsertsAndRemovesKnot0_5) {  // NOLINT
  nurbs_2d_after_->InsertKnot(ParamCoord(0.5), 0);
  ASSERT_THAT(nurbs_2d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_2d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  nurbs_2d_after_->RemoveKnot(ParamCoord(0.5), 0, 1e-15);
  ASSERT_THAT(nurbs_2d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_2d_before_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(nurbs_2d_after_->GetNumberOfControlPoints(), nurbs_2d_before_->GetNumberOfControlPoints());
  for (int point = 0; point < nurbs_2d_before_->GetNumberOfControlPoints(); ++point) {
    for (int dimension = 0; dimension < nurbs_2d_before_->GetDimension(); ++dimension) {
      ASSERT_THAT(nurbs_2d_after_->GetControlPoint({point}, dimension),
                  DoubleNear(nurbs_2d_before_->GetControlPoint({point}, dimension), 1e-12));
    }
    ASSERT_THAT(nurbs_2d_after_->GetWeight({point}), DoubleNear(nurbs_2d_before_->GetWeight({point}), 1e-12));
  }
  std::array<ParamCoord, 2> param_coord{};
  for (int i = 0; i < 10; ++i) {
    for (int coord = 0; coord < 2; ++coord) {
      param_coord[coord] = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    }
    for (int dimension = 0; dimension < nurbs_2d_before_->GetDimension(); ++dimension) {
      ASSERT_THAT(nurbs_2d_after_->Evaluate(param_coord, {dimension})[0],
                  DoubleNear(nurbs_2d_before_->Evaluate(param_coord, {dimension})[0], 1e-12));
    }
  }
}

class Random3DNURBSForKnotRemoval : public Test {  // NOLINT
 public:
  Random3DNURBSForKnotRemoval() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0.0}, ParamCoord{1.0}};
    spl::RandomNURBSGenerator<3> spline_generator(limits, 5, 3);
    spl::NURBS<3> nurbs(spline_generator);
    nurbs_3d_before_ = std::make_shared<spl::NURBS<3>>(nurbs);
    spl::NURBS<3> nurbs_after(nurbs);
    nurbs_3d_after_ = std::make_shared<spl::NURBS<3>>(nurbs_after);
  }

 protected:
  std::shared_ptr<spl::NURBS<3>> nurbs_3d_before_;
  std::shared_ptr<spl::NURBS<3>> nurbs_3d_after_;
};

TEST_F(Random3DNURBSForKnotRemoval, InsertsAndRemovesKnot0_5) {  // NOLINT
  nurbs_3d_after_->InsertKnot(ParamCoord(0.5), 0);
  ASSERT_THAT(nurbs_3d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_3d_before_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  nurbs_3d_after_->RemoveKnot(ParamCoord(0.5), 0, 1e-15);
  ASSERT_THAT(nurbs_3d_after_->GetKnotVector(0)->GetNumberOfKnots(),
              nurbs_3d_before_->GetKnotVector(0)->GetNumberOfKnots());
  ASSERT_THAT(nurbs_3d_after_->GetNumberOfControlPoints(), nurbs_3d_before_->GetNumberOfControlPoints());
  for (int point = 0; point < nurbs_3d_before_->GetNumberOfControlPoints(); ++point) {
    for (int dimension = 0; dimension < nurbs_3d_before_->GetDimension(); ++dimension) {
      ASSERT_THAT(nurbs_3d_after_->GetControlPoint({point}, dimension),
                  DoubleNear(nurbs_3d_before_->GetControlPoint({point}, dimension), 1e-12));
    }
    ASSERT_THAT(nurbs_3d_after_->GetWeight({point}), DoubleNear(nurbs_3d_before_->GetWeight({point}), 1e-12));
  }
  std::array<ParamCoord, 3> param_coord{};
  for (int i = 0; i < 10; ++i) {
    for (int coord = 0; coord < 3; ++coord) {
      param_coord[coord] = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    }
    for (int dimension = 0; dimension < nurbs_3d_before_->GetDimension(); ++dimension) {
      ASSERT_THAT(nurbs_3d_after_->Evaluate(param_coord, {dimension})[0],
                  DoubleNear(nurbs_3d_before_->Evaluate(param_coord, {dimension})[0], 1e-12));
    }
  }
}
