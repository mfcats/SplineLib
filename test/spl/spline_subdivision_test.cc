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
using testing::DoubleEq;
using testing::DoubleNear;

class BSplineFig5_9 : public Test {  // NOLINT
 public:
  BSplineFig5_9() {
    std::array<Degree, 2> degree = {Degree{3}, Degree{2}};
    std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector = {
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                           ParamCoord{1}, ParamCoord{1}, ParamCoord{1},
                                                           ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5},
                                                           ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({5.0, 0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 3.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 0.0, 3.0})),

        baf::ControlPoint(std::vector<double>({5.0, 2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.0, 1.5, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 3.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 1.0, 3.0})),

        baf::ControlPoint(std::vector<double>({5.0, 2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 2.5, 2.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 2.5, 2.0})),

        baf::ControlPoint(std::vector<double>({5.0, 5.0, 0.0})),
        baf::ControlPoint(std::vector<double>({2.0, 4.5, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 4.0, 2.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 4.0, 2.0}))
    };
    bspline_2d_ = std::make_shared<spl::BSpline<2>>(knot_vector, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> bspline_2d_;
};

TEST_F(BSplineFig5_9, IsSubdividedAtKnot0_3InFirstDirection) {  // NOLINT
  auto splines = bspline_2d_->SudivideSpline(ParamCoord{0.7}, 0);
  auto spline1 = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(splines[0]);
  auto spline2 = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(splines[1]);
  ASSERT_THAT(spline1->GetKnotVector(0)->GetNumberOfKnots(), 8);
  ASSERT_THAT(spline2->GetKnotVector(0)->GetNumberOfKnots(), 8);
  ASSERT_THAT(spline1->GetKnotVector(0)->GetKnot(0).get(), DoubleEq(0));
  ASSERT_THAT(spline1->GetKnotVector(0)->GetKnot(3).get(), DoubleEq(0));
  ASSERT_THAT(spline1->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(0.7));
  ASSERT_THAT(spline1->GetKnotVector(0)->GetKnot(7).get(), DoubleEq(0.7));
  ASSERT_THAT(spline2->GetKnotVector(0)->GetKnot(0).get(), DoubleEq(0.7));
  ASSERT_THAT(spline2->GetKnotVector(0)->GetKnot(3).get(), DoubleEq(0.7));
  ASSERT_THAT(spline2->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(1));
  ASSERT_THAT(spline2->GetKnotVector(0)->GetKnot(7).get(), DoubleEq(1));
  for (int i = 0; i <= 100; ++i) {
    for (int j = 0; j <= 100; ++j) {
      std::array<ParamCoord, 2> param_coord{ParamCoord(i / 100.0), ParamCoord(j / 100.0)};
      if (i / 100.0 <= 0.7) {
        ASSERT_THAT(spline1->Evaluate(param_coord, {0})[0],
                    DoubleNear(bspline_2d_->Evaluate(param_coord, {0})[0], 0.000001));
      }
      if (i / 100.0 >= 0.7) {
        ASSERT_THAT(spline2->Evaluate(param_coord, {0})[0],
                    DoubleNear(bspline_2d_->Evaluate(param_coord, {0})[0], 0.000001));
      }
    }
  }
}

class Random1DBSplineToSplit : public Test {  // NOLINT
 public:
  Random1DBSplineToSplit() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0}, ParamCoord{1}};
    spl::RandomBSplineGenerator<1> spline_generator(limits, 10, 2);
    spl::BSpline<1> b_spline(*spline_generator.GetParameterSpace(), *spline_generator.GetPhysicalSpace());
    b_spline_1d_ = std::make_shared<spl::BSpline<1>>(b_spline);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> b_spline_1d_;
};

TEST_F(Random1DBSplineToSplit, IsSubdividedAtKnot0_25) {  // NOLINT
  auto splines = b_spline_1d_->SudivideSpline(ParamCoord{0.25}, 0);
  auto spline1 = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(splines[0]);
  auto spline2 = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(splines[1]);
  double eval_points = 100;
  for (int i = 0; i <= eval_points; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / eval_points)};
    if (i / eval_points <= 0.25) {
      ASSERT_THAT(spline1->Evaluate(param_coord, {0})[0],
                  DoubleNear(b_spline_1d_->Evaluate(param_coord, {0})[0], 0.000001));
    }
    if (i / eval_points >= 0.25) {
      ASSERT_THAT(spline2->Evaluate(param_coord, {0})[0],
                  DoubleNear(b_spline_1d_->Evaluate(param_coord, {0})[0], 0.000001));
    }
  }
}

class Random1DNURBSToSplit : public Test {  // NOLINT
 public:
  Random1DNURBSToSplit() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0}, ParamCoord{1}};
    spl::RandomNURBSGenerator<1> spline_generator(limits, 10, 3);
    spl::NURBS<1> nurbs(spline_generator.GetParameterSpace(), spline_generator.GetWeightedPhysicalSpace());
    nurbs_1d_ = std::make_shared<spl::NURBS<1>>(nurbs);
  }

 protected:
  std::shared_ptr<spl::NURBS<1>> nurbs_1d_;
};

TEST_F(Random1DNURBSToSplit, IsSubdividedAtKnot0_99) {  // NOLINT
  auto splines = nurbs_1d_->SudivideSpline(ParamCoord{0.99}, 0);
  auto spline1 = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(splines[0]);
  auto spline2 = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(splines[1]);
  double eval_points = 100;
  for (int i = 0; i <= eval_points; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / eval_points)};
    if (i / eval_points <= 0.99) {
      ASSERT_THAT(spline1->Evaluate(param_coord, {0})[0],
                  DoubleNear(nurbs_1d_->Evaluate(param_coord, {0})[0], 0.000001));
    }
    if (i / eval_points >= 0.99) {
      ASSERT_THAT(spline2->Evaluate(param_coord, {0})[0],
                  DoubleNear(nurbs_1d_->Evaluate(param_coord, {0})[0], 0.000001));
    }
  }
}

class Random2DBSplineToSplit : public Test {  // NOLINT
 public:
  Random2DBSplineToSplit() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0}, ParamCoord{1}};
    spl::RandomBSplineGenerator<2> spline_generator(limits, 10, 3);
    spl::BSpline<2> b_spline(*spline_generator.GetParameterSpace(), *spline_generator.GetPhysicalSpace());
    b_spline_2d_ = std::make_shared<spl::BSpline<2>>(b_spline);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> b_spline_2d_;
};

TEST_F(Random2DBSplineToSplit, IsSubdividedAtKnot0_1InFirstDirection) {  // NOLINT
  auto splines = b_spline_2d_->SudivideSpline(ParamCoord{0.1}, 0);
  auto spline1 = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(splines[0]);
  auto spline2 = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(splines[1]);
  double eval_points = 25;
  for (int i = 0; i <= eval_points; ++i) {
    for (int j = 0; j <= eval_points; ++j) {
      std::array<ParamCoord, 2> param_coord{ParamCoord(i / eval_points), ParamCoord(j / eval_points)};
      if (i / eval_points <= 0.1) {
        ASSERT_THAT(spline1->Evaluate(param_coord, {0})[0],
                    DoubleNear(b_spline_2d_->Evaluate(param_coord, {0})[0], 0.000001));
      }
      if (i / eval_points >= 0.1) {
        ASSERT_THAT(spline2->Evaluate(param_coord, {0})[0],
                    DoubleNear(b_spline_2d_->Evaluate(param_coord, {0})[0], 0.000001));
      }
    }
  }
}

class Random2DNURBSToSplit : public Test {  // NOLINT
 public:
  Random2DNURBSToSplit() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0}, ParamCoord{1}};
    spl::RandomNURBSGenerator<2> spline_generator(limits, 10, 3);
    spl::NURBS<2> nurbs(spline_generator.GetParameterSpace(), spline_generator.GetWeightedPhysicalSpace());
    nurbs_2d_ = std::make_shared<spl::NURBS<2>>(nurbs);
  }

 protected:
  std::shared_ptr<spl::NURBS<2>> nurbs_2d_;
};

TEST_F(Random2DNURBSToSplit, IsSubdividedAtKnot0_5InSecondDirection) {  // NOLINT
  auto splines = nurbs_2d_->SudivideSpline(ParamCoord{0.5}, 1);
  auto spline1 = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(splines[0]);
  auto spline2 = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(splines[1]);
  double eval_points = 25;
  for (int i = 0; i <= eval_points; ++i) {
    for (int j = 0; j <= eval_points; ++j) {
      std::array<ParamCoord, 2> param_coord{ParamCoord(i / eval_points), ParamCoord(j / eval_points)};
      if (j / eval_points <= 0.5) {
        ASSERT_THAT(spline1->Evaluate(param_coord, {0})[0],
                    DoubleNear(nurbs_2d_->Evaluate(param_coord, {0})[0], 0.000001));
      }
      if (j / eval_points >= 0.5) {
        ASSERT_THAT(spline2->Evaluate(param_coord, {0})[0],
                    DoubleNear(nurbs_2d_->Evaluate(param_coord, {0})[0], 0.000001));
      }
    }
  }
}

class Random3DBSplineToSplit : public Test {  // NOLINT
 public:
  Random3DBSplineToSplit() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0}, ParamCoord{1}};
    spl::RandomBSplineGenerator<3> spline_generator(limits, 10, 3);
    spl::BSpline<3> b_spline(*spline_generator.GetParameterSpace(), *spline_generator.GetPhysicalSpace());
    b_spline_3d_ = std::make_shared<spl::BSpline<3>>(b_spline);
  }

 protected:
  std::shared_ptr<spl::BSpline<3>> b_spline_3d_;
};

TEST_F(Random3DBSplineToSplit, IsSubdividedAtKnot0_4InSecondDirection) {  // NOLINT
  auto splines = b_spline_3d_->SudivideSpline(ParamCoord{0.4}, 1);
  auto spline1 = std::any_cast<std::shared_ptr<spl::BSpline<3>>>(splines[0]);
  auto spline2 = std::any_cast<std::shared_ptr<spl::BSpline<3>>>(splines[1]);
  double eval_points = 15;
  for (int i = 0; i <= eval_points; ++i) {
    for (int j = 0; j <= eval_points; ++j) {
      for (int k = 0; k <= eval_points; ++k) {
        std::array<ParamCoord, 3>
            param_coord{ParamCoord(i / eval_points), ParamCoord(j / eval_points), ParamCoord(k / eval_points)};
        if (j / eval_points <= 0.4) {
          ASSERT_THAT(spline1->Evaluate(param_coord, {0})[0],
                      DoubleNear(b_spline_3d_->Evaluate(param_coord, {0})[0], 0.000001));
        }
        if (j / eval_points >= 0.4) {
          ASSERT_THAT(spline2->Evaluate(param_coord, {0})[0],
                      DoubleNear(b_spline_3d_->Evaluate(param_coord, {0})[0], 0.000001));
        }
      }
    }
  }
}

class Random3DNURBSToSplit : public Test {  // NOLINT
 public:
  Random3DNURBSToSplit() {
    std::array<ParamCoord, 2> limits = {ParamCoord{0}, ParamCoord{1}};
    spl::RandomNURBSGenerator<3> spline_generator(limits, 10, 3);
    spl::NURBS<3> nurbs(spline_generator.GetParameterSpace(), spline_generator.GetWeightedPhysicalSpace());
    nurbs_3d_ = std::make_shared<spl::NURBS<3>>(nurbs);
  }

 protected:
  std::shared_ptr<spl::NURBS<3>> nurbs_3d_;
};

TEST_F(Random3DNURBSToSplit, IsSubdividedAtKnot0_9InThirdDirection) {  // NOLINT
  auto splines = nurbs_3d_->SudivideSpline(ParamCoord{0.9}, 2);
  auto spline1 = std::any_cast<std::shared_ptr<spl::NURBS<3>>>(splines[0]);
  auto spline2 = std::any_cast<std::shared_ptr<spl::NURBS<3>>>(splines[1]);
  double s = 5;
  for (int i = 0; i <= s; ++i) {
    for (int j = 0; j <= s; ++j) {
      for (int k = 0; k <= s; ++k) {
        std::array<ParamCoord, 3> param_coord{ParamCoord(i / s), ParamCoord(j / s), ParamCoord(k / s)};
        if (k / s <= 0.9) {
          ASSERT_THAT(spline1->Evaluate(param_coord, {0})[0],
                      DoubleNear(nurbs_3d_->Evaluate(param_coord, {0})[0], 0.000001));
        }
        if (k / s >= 0.9) {
          ASSERT_THAT(spline2->Evaluate(param_coord, {0})[0],
                      DoubleNear(nurbs_3d_->Evaluate(param_coord, {0})[0], 0.000001));
        }
      }
    }
  }
}
