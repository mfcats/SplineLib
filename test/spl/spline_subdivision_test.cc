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

TEST_F(BSplineFig5_9, IsSubdividedAtKnot0_7InFirstDirection) {  // NOLINT
  auto splines = bspline_2d_->SudivideSpline(ParamCoord{0.7}, 0);
  ASSERT_THAT(splines[0]->GetKnotVector(0)->GetNumberOfKnots(), 8);
  ASSERT_THAT(splines[1]->GetKnotVector(0)->GetNumberOfKnots(), 8);
  ASSERT_THAT(splines[0]->GetKnotVector(0)->GetKnot(0).get(), DoubleEq(0));
  ASSERT_THAT(splines[0]->GetKnotVector(0)->GetKnot(3).get(), DoubleEq(0));
  ASSERT_THAT(splines[0]->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(0.7));
  ASSERT_THAT(splines[0]->GetKnotVector(0)->GetKnot(7).get(), DoubleEq(0.7));
  ASSERT_THAT(splines[1]->GetKnotVector(0)->GetKnot(0).get(), DoubleEq(0.7));
  ASSERT_THAT(splines[1]->GetKnotVector(0)->GetKnot(3).get(), DoubleEq(0.7));
  ASSERT_THAT(splines[1]->GetKnotVector(0)->GetKnot(4).get(), DoubleEq(1));
  ASSERT_THAT(splines[1]->GetKnotVector(0)->GetKnot(7).get(), DoubleEq(1));
  double steps = 100;
  for (int i = 0; i <= steps; ++i) {
    ParamCoord coord2 = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    std::array<ParamCoord, 2> param_coord{ParamCoord(i / steps), coord2};
    int spline_number = i / steps >= 0.7 ? 1 : 0;
    std::vector<double> evaluated_splitted_spline = splines[spline_number]->Evaluate(param_coord, {0, 1});
    std::vector<double> original_spline = bspline_2d_->Evaluate(param_coord, {0, 1});
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(evaluated_splitted_spline[j], DoubleNear(original_spline[j], 0.000001));
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
  double steps = 100;
  for (int i = 0; i <= steps; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / steps)};
    int spline_number = i / steps >= 0.25 ? 1 : 0;
    double evaluated_splitted_spline = splines[spline_number]->Evaluate(param_coord, {0})[0];
    double original_spline = b_spline_1d_->Evaluate(param_coord, {0})[0];
    ASSERT_THAT(evaluated_splitted_spline, DoubleNear(original_spline, 0.000001));
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
  double steps = 100;
  for (int i = 0; i <= steps; ++i) {
    std::array<ParamCoord, 1> param_coord{ParamCoord(i / steps)};
    int spline_number = i / steps >= 0.99 ? 1 : 0;
    double evaluated_splitted_spline = splines[spline_number]->Evaluate(param_coord, {0})[0];
    double original_spline = nurbs_1d_->Evaluate(param_coord, {0})[0];
    ASSERT_THAT(evaluated_splitted_spline, DoubleNear(original_spline, 0.000001));
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
  double steps = 50;
  for (int i = 0; i <= steps; ++i) {
    ParamCoord coord2 = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    std::array<ParamCoord, 2> param_coord{ParamCoord(i / steps), coord2};
    int spline_number = i / steps >= 0.1 ? 1 : 0;
    std::vector<double> evaluated_splitted_spline = splines[spline_number]->Evaluate(param_coord, {0, 1});
    std::vector<double> original_spline = b_spline_2d_->Evaluate(param_coord, {0, 1});
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(evaluated_splitted_spline[j], DoubleNear(original_spline[j], 0.000001));
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
  double steps = 50;
  for (int i = 0; i <= steps; ++i) {
    ParamCoord coord1 = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    std::array<ParamCoord, 2> param_coord{coord1, ParamCoord(i / steps)};
    int spline_number = i / steps >= 0.5 ? 1 : 0;
    std::vector<double> evaluated_splitted_spline = splines[spline_number]->Evaluate(param_coord, {0, 1});
    std::vector<double> original_spline = nurbs_2d_->Evaluate(param_coord, {0, 1});
    for (int j = 0; j < 2; ++j) {
      ASSERT_THAT(evaluated_splitted_spline[j], DoubleNear(original_spline[j], 0.000001));
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
  double steps = 25;
  for (int i = 0; i <= steps; ++i) {
    ParamCoord coord1 = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    ParamCoord coord3 = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    std::array<ParamCoord, 3> param_coord{coord1, ParamCoord(i / steps), coord3};
    int spline_number = i / steps >= 0.4 ? 1 : 0;
    std::vector<double> evaluated_splitted_spline = splines[spline_number]->Evaluate(param_coord, {0, 1, 2});
    std::vector<double> original_spline = b_spline_3d_->Evaluate(param_coord, {0, 1, 2});
    for (int j = 0; j < 3; ++j) {
      ASSERT_THAT(evaluated_splitted_spline[j], DoubleNear(original_spline[j], 0.000001));
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
  double steps = 5;
  for (int i = 0; i <= steps; ++i) {
    ParamCoord coord1 = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    ParamCoord coord2 = ParamCoord{util::Random::GetUniformRandom<double>(0.0, 1.0)};
    std::array<ParamCoord, 3> param_coord{coord1, coord2, ParamCoord(i / steps)};
    int spline_number = i / steps >= 0.9 ? 1 : 0;
    std::vector<double> evaluated_splitted_spline = splines[spline_number]->Evaluate(param_coord, {0, 1, 2});
    std::vector<double> original_spline = nurbs_3d_->Evaluate(param_coord, {0, 1, 2});
    for (int j = 0; j < 3; ++j) {
      ASSERT_THAT(evaluated_splitted_spline[j], DoubleNear(original_spline[j], 0.000001));
    }
  }
}
