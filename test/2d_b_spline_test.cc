/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "b_spline.h"

#include "gmock/gmock.h"

using testing::Test;

class A2DBSpline : public Test {
 public:
  A2DBSpline() {
    std::array<baf::KnotVector, 2> knot_vector =
        {baf::KnotVector({std::vector<ParamCoord>({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                   ParamCoord{1}, ParamCoord{1}})}),
         baf::KnotVector({std::vector<ParamCoord>({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                   ParamCoord{1}, ParamCoord{1}})})};
    std::array<Degree, 2> degree = {Degree{2}, Degree{2}};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({-1.0, -1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, -1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, -1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0}))
    };
    knot_vector_ptr = std::make_shared<std::array<baf::KnotVector, 2>>(knot_vector);
    b_spline = std::make_unique<spl::BSpline<2>>(knot_vector_ptr, degree, control_points);
  }

 protected:
  std::unique_ptr<spl::BSpline<2>> b_spline;
  std::shared_ptr<std::array<baf::KnotVector, 2>> knot_vector_ptr;
};

TEST_F(A2DBSpline, Corner) { // NOLINT
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {0})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {1})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {2})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, EdgeDim0) { // NOLINT
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.0}, ParamCoord{0.33333}}, {0})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.0}, ParamCoord{0.33333}}, {1})[0], -0.33333, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.0}, ParamCoord{0.33333}}, {2})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, EdgeDim1) { // NOLINT
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.33333}, ParamCoord{0.0}}, {0})[0], -0.33333, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.33333}, ParamCoord{0.0}}, {1})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.33333}, ParamCoord{0.0}}, {2})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, Center) { // NOLINT
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.5}, ParamCoord{0.5}}, {0})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.5}, ParamCoord{0.5}}, {1})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.5}, ParamCoord{0.5}}, {2})[0], 0.25, 0.00005);
}

TEST_F(A2DBSpline, Random) { // NOLINT
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.75}, ParamCoord{0.25}}, {0})[0], 0.5, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.75}, ParamCoord{0.25}}, {1})[0], -0.5, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.75}, ParamCoord{0.25}}, {2})[0], 0.14063, 0.00005);
}

TEST_F(A2DBSpline, CornerDer10) { // NOLINT
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.0}}, {0}, {1, 0})[0], 2.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.0}}, {1}, {1, 0})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.0}}, {2}, {1, 0})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, CornerDer01) { // NOLINT
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.0}}, {0}, {0, 1})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.0}}, {1}, {0, 1})[0], 2.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.0}}, {2}, {0, 1})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, EdgeDim0Der10) { // NOLINT
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.33333}}, {0}, {1, 0})[0], 2.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.33333}}, {1}, {1, 0})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.33333}}, {2}, {1, 0})[0], 0.888889, 0.00005);
}

TEST_F(A2DBSpline, EdgeDim0Der01) { // NOLINT
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.33333}}, {0}, {0, 1})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.33333}}, {1}, {0, 1})[0], 2.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.33333}}, {2}, {0, 1})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, CenterDer10) { // NOLINT
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}}, {0}, {1, 0})[0], 2.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}}, {1}, {1, 0})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}}, {2}, {1, 0})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, RandomDer10) { // NOLINT
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {0}, {1, 0})[0], 2.000, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {1}, {1, 0})[0], 0.000, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {2}, {1, 0})[0], -0.375, 0.00005);
}

TEST_F(A2DBSpline, RandomDer01) { // NOLINT
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {0}, {0, 1})[0], 0.000, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {1}, {0, 1})[0], 2.000, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {2}, {0, 1})[0], 0.375, 0.00005);
}

TEST_F(A2DBSpline, RandomDer12) { // NOLINT
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {0}, {1, 2})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {1}, {1, 2})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {2}, {1, 2})[0], 4.0, 0.00005);
}
