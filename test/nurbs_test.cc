/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <numeric_settings.h>
#include "gmock/gmock.h"

#include "nurbs.h"

using testing::Test;
using testing::DoubleNear;

class NurbsEx4_1 : public Test {
 public:
  NurbsEx4_1() {
    std::array<baf::KnotVector, 1> knot_vector = {baf::KnotVector({0, 0, 0, 1, 2, 3, 3, 3})};
    std::array<int, 1> degree = {2};
    std::vector<double> weights = {1, 4, 1, 1, 1};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({5.0, -1.0}))
    };
    nurbs = std::make_unique<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::unique_ptr<spl::NURBS<1>> nurbs;
};

TEST_F(NurbsEx4_1, Returns1_4For1AndDim0) {
  ASSERT_THAT(nurbs->Evaluate({1.0}, {0})[0], DoubleNear(1.4, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(NurbsEx4_1, Returns1_2For1AndDim1) {
  ASSERT_THAT(nurbs->Evaluate({1.0}, {1})[0], DoubleNear(1.2, util::NumericSettings<double>::kEpsilon()));
}

class ANurbs : public Test {
 public:
  ANurbs() {
    std::array<baf::KnotVector, 1>
        knot_vector = {baf::KnotVector({0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 0.95, 1.0, 1.0, 1.0})};
    std::array<int, 1> degree = {2};
    std::vector<double> weights = {1.0, 0.9, 0.7, 0.5, 0.8, 1.2, 2.0};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.5, 3.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.5, 5.5, 4.0})),
        baf::ControlPoint(std::vector<double>({4.5, 5.5, 0.1})),
        baf::ControlPoint(std::vector<double>({3.0, 1.5, 2.0})),
        baf::ControlPoint(std::vector<double>({7.5, 1.5, 3.5})),
        baf::ControlPoint(std::vector<double>({6.0, 4.0, 5.3})),
        baf::ControlPoint(std::vector<double>({8.5, 4.5, 0.0}))
    };
    nurbs = std::make_unique<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::unique_ptr<spl::NURBS<1>> nurbs;
};

TEST_F(ANurbs, ReturnsCorrectCurvePointForFirstKnot) {
  ASSERT_THAT(nurbs->Evaluate({0.0}, {0})[0], DoubleNear(0.5, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({0.0}, {1})[0], DoubleNear(3.0, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({0.0}, {2})[0], DoubleNear(1.0, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(ANurbs, ReturnsCorrectCurvePointForInnerKnot) {
  ASSERT_THAT(nurbs->Evaluate({0.25}, {0})[0], DoubleNear(2.8125, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({0.25}, {1})[0], DoubleNear(5.5, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({0.25}, {2})[0], DoubleNear(2.29375, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(ANurbs, ReturnsCorrectCurvePointForValueBetweenTwoKnots) {
  ASSERT_THAT(nurbs->Evaluate({1.0 / 3.0}, {0})[0], DoubleNear(3.625, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({1.0 / 3.0}, {1})[0], DoubleNear(5.34848, 0.000005));
  ASSERT_THAT(nurbs->Evaluate({1.0 / 3.0}, {2})[0], DoubleNear(1.23561, 0.000005));
}

TEST_F(ANurbs, ReturnsCorrectCurvePointForLastKnot) {
  ASSERT_THAT(nurbs->Evaluate({1.0}, {0})[0], DoubleNear(8.5, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({1.0}, {1})[0], DoubleNear(4.5, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({1.0}, {2})[0], DoubleNear(0.0, util::NumericSettings<double>::kEpsilon()));
}

class NurbsDerivativeEx4_2 : public Test {
 public:
  NurbsDerivativeEx4_2() {
    std::array<baf::KnotVector, 1> knot_vector = {baf::KnotVector({0, 0, 0, 1, 1, 1})};
    std::array<int, 1> degree = {2};
    std::vector<double> weights = {1, 1, 2};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0}))
    };
    nurbs = std::make_unique<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::unique_ptr<spl::NURBS<1>> nurbs;
};

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForFirstDerivativeAtFirstKnot) {
  ASSERT_THAT(nurbs->EvaluateDerivative({0.0}, {0}, {1})[0], 0.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({0.0}, {1}, {1})[0], 2.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForFirstDerivativeAtValueBetweenKnots) {
  ASSERT_THAT(nurbs->EvaluateDerivative({0.5}, {0}, {1})[0], -1.28);
  ASSERT_THAT(nurbs->EvaluateDerivative({0.5}, {1}, {1})[0], 0.96);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForFirstDerivativeAtLastKnot) {
  ASSERT_THAT(nurbs->EvaluateDerivative({1.0}, {0}, {1})[0], -1.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({1.0}, {1}, {1})[0], 0.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForSecondDerivativeAtFirstKnot) {
  ASSERT_THAT(nurbs->EvaluateDerivative({0.0}, {0}, {2})[0], -4.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({0.0}, {1}, {2})[0], 0.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForSecondDerivativeAtValueBetweenKnots) {
  ASSERT_THAT(nurbs->EvaluateDerivative({0.5}, {0}, {2})[0],
              DoubleNear(-0.512, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->EvaluateDerivative({0.5}, {1}, {2})[0], -2.816);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForSecondDerivativeAtLastKnot) {
  ASSERT_THAT(nurbs->EvaluateDerivative({1.0}, {0}, {2})[0], 1.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({1.0}, {1}, {2})[0], -1.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForThirdDerivativeAtFirstKnot) {
  ASSERT_THAT(nurbs->EvaluateDerivative({0.0}, {0}, {3})[0], 0.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({0.0}, {1}, {3})[0], -12.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForThirdDerivativeAtValueBetweenKnots) {
  ASSERT_THAT(nurbs->EvaluateDerivative({0.5}, {0}, {3})[0],
              DoubleNear(7.3728, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->EvaluateDerivative({0.5}, {1}, {3})[0],
              DoubleNear(2.1504, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForThirdDerivativeAtLastKnot) {
  ASSERT_THAT(nurbs->EvaluateDerivative({1.0}, {0}, {3})[0], 0.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({1.0}, {1}, {3})[0], 3.0);
}
