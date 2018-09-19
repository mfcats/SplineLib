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
#include "nurbs_generator.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;
using ::testing::NiceMock;

class MockParameterSpace : public spl::ParameterSpace<1> {
 public:
  MOCK_CONST_METHOD1(GetDegree, Degree(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 1>, std::array<ParamCoord, 1>));
  MOCK_CONST_METHOD3(GetBasisFunctionDerivatives,
                     double(std::array<int, 1>, std::array<ParamCoord, 1>, std::array<int, 1>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 1>(std::array<ParamCoord, 1>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParamCoord, 1>));
};

class MockWeightedPhysicalSpace : public spl::WeightedPhysicalSpace<1> {
 public:
  MOCK_CONST_METHOD1(GetControlPoint, baf::ControlPoint(std::array<int, 1>));
};

class NurbsEx4_1 : public Test {
 public:
  NurbsEx4_1() : 
      parameter_space(std::make_shared<NiceMock<MockParameterSpace>>()),
      weihted_physical_space(std::make_shared<NiceMock<MockWeightedPhysicalSpace>>()) {
    spl::NURBSGenerator<1> nurbs_generator(weighted_physical_space, parameter_space);
    nurbs = std::make_unique<spl::NURBS<1>>(nurbs_generator);
    /*
    std::array<baf::KnotVector, 1> knot_vector =
    {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2}, ParamCoord{3},
                          ParamCoord{3}, ParamCoord{3}})};
    std::array<Degree, 1> degree = {Degree{2}};
    std::vector<double> weights = {1, 4, 1, 1, 1};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({5.0, -1.0}))
    };
    std::array<std::shared_ptr<baf::KnotVector>, 1>
        knot_vector_ptr = {std::make_shared<baf::KnotVector>(knot_vector[0])};
    nurbs = std::make_unique<spl::NURBS<1>>(knot_vector_ptr, degree, control_points, weights);
    */
  }

 protected:
  std::unique_ptr<spl::NURBS<1>> nurbs;
  std::shared_ptr<NiceMock<MockParameterSpace>> parameter_space;
  std::shared_ptr<NiceMock<MockWeightedPhysicalSpace>> weighted_physical_space;
};

TEST_F(NurbsEx4_1, Returns1_4For1AndDim0) { // NOLINT
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{1.0}}, {0})[0], DoubleNear(1.4, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(NurbsEx4_1, Returns1_2For1AndDim1) { // NOLINT
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{1.0}}, {1})[0], DoubleNear(1.2, util::NumericSettings<double>::kEpsilon()));
}

class ANurbs : public Test {
 public:
  ANurbs() : degree_{Degree{2}} {
    std::array<baf::KnotVector, 1>
        knot_vector =
        {baf::KnotVector({ParamCoord{0.0}, ParamCoord{0.0}, ParamCoord{0.0}, ParamCoord{0.25}, ParamCoord{0.5},
                          ParamCoord{0.75}, ParamCoord{0.95}, ParamCoord{1.0}, ParamCoord{1.0}, ParamCoord{1.0}})};
    weights_ = {1.0, 0.9, 0.7, 0.5, 0.8, 1.2, 2.0};
    control_points_ = {
        baf::ControlPoint(std::vector<double>({0.5, 3.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.5, 5.5, 4.0})),
        baf::ControlPoint(std::vector<double>({4.5, 5.5, 0.1})),
        baf::ControlPoint(std::vector<double>({3.0, 1.5, 2.0})),
        baf::ControlPoint(std::vector<double>({7.5, 1.5, 3.5})),
        baf::ControlPoint(std::vector<double>({6.0, 4.0, 5.3})),
        baf::ControlPoint(std::vector<double>({8.5, 4.5, 0.0}))
    };
    knot_vector_[0] = std::make_shared<baf::KnotVector>(knot_vector[0]);
    nurbs = std::make_unique<spl::NURBS<1>>(knot_vector_, degree_, control_points_, weights_);
  }

 protected:
  std::unique_ptr<spl::NURBS<1>> nurbs;
  std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_;
  std::array<Degree, 1> degree_;
  std::vector<double> weights_;
  std::vector<baf::ControlPoint> control_points_;
};

TEST_F(ANurbs, ReturnsCorrectCurvePointForFirstKnot) { // NOLINT
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{0.0}}, {0})[0], DoubleNear(0.5, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{0.0}}, {1})[0], DoubleNear(3.0, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{0.0}}, {2})[0], DoubleNear(1.0, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(ANurbs, ReturnsCorrectCurvePointForInnerKnot) { // NOLINT
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{0.25}}, {0})[0],
              DoubleNear(2.8125, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{0.25}}, {1})[0], DoubleNear(5.5, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{0.25}}, {2})[0],
              DoubleNear(2.29375, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(ANurbs, ReturnsCorrectCurvePointForValueBetweenTwoKnots) { // NOLINT
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{1.0 / 3.0}}, {0})[0],
              DoubleNear(3.625, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{1.0 / 3.0}}, {1})[0], DoubleNear(5.34848, 0.000005));
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{1.0 / 3.0}}, {2})[0], DoubleNear(1.23561, 0.000005));
}

TEST_F(ANurbs, ReturnsCorrectCurvePointForLastKnot) { // NOLINT
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{1.0}}, {0})[0], DoubleNear(8.5, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{1.0}}, {1})[0], DoubleNear(4.5, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{1.0}}, {2})[0], DoubleNear(0.0, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(ANurbs, CanBeConstructedWithAPhysicalAndAParameterSpace) { // NOLINT
  std::shared_ptr<spl::ParameterSpace<1>>
      parameter_space = std::make_shared<spl::ParameterSpace<1>>(spl::ParameterSpace<1>(knot_vector_, {degree_[0]}));
  std::shared_ptr<spl::WeightedPhysicalSpace<1>> physicalSpace =
      std::make_shared<spl::WeightedPhysicalSpace<1>>(spl::WeightedPhysicalSpace<1>(control_points_, weights_, {7}));
  nurbs = std::make_unique<spl::NURBS<1>>(parameter_space, physicalSpace);
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{0.0}}, {0})[0], DoubleEq(0.5));
}

TEST_F(ANurbs, ThrowsExceptionForDifferentControlPointDimensions) { // NOLINT
  control_points_.emplace_back(std::vector<double>({0.0}));
  ASSERT_THROW(std::make_unique<spl::NURBS<1>>(knot_vector_, degree_, control_points_, weights_), std::runtime_error);
}

TEST_F(ANurbs, ThrowsExceptionForEvaluationAt1_2) { // NOLINT
  ASSERT_THROW(nurbs->Evaluate({ParamCoord{1.2}}, {0}), std::runtime_error);
}

TEST_F(ANurbs, ThrowsExceptionForEvaluationAtMinus0_1) { // NOLINT
  ASSERT_THROW(nurbs->Evaluate({ParamCoord{-0.1}}, {0}), std::runtime_error);
}

TEST_F(ANurbs, ThrowsExceptionForDerivativeEvaluationAt1_05) { // NOLINT
  ASSERT_THROW(nurbs->EvaluateDerivative({ParamCoord{1.05}}, {0}, {1}), std::runtime_error);
}

TEST_F(ANurbs, ThrowsExceptionForDerivativeEvaluationAtMinus1_0) { // NOLINT
  ASSERT_THROW(nurbs->EvaluateDerivative({ParamCoord{-0.05}}, {0}, {1}), std::runtime_error);
}

class NurbsDerivativeEx4_2 : public Test {
 public:
  NurbsDerivativeEx4_2() {
    std::array<baf::KnotVector, 1> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})};
    std::array<Degree, 1> degree = {Degree{2}};
    std::vector<double> weights = {1, 1, 2};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0}))
    };
    std::array<std::shared_ptr<baf::KnotVector>, 1>
        knot_vector_ptr = {std::make_shared<baf::KnotVector>(knot_vector[0])};
    nurbs = std::make_unique<spl::NURBS<1>>(knot_vector_ptr, degree, control_points, weights);
  }

 protected:
  std::unique_ptr<spl::NURBS<1>> nurbs;
};

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForFirstDerivativeAtFirstKnot) { // NOLINT
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{0.0}}, {0}, {1})[0], 0.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{0.0}}, {1}, {1})[0], 2.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForFirstDerivativeAtValueBetweenKnots) { // NOLINT
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{0.5}}, {0}, {1})[0], -1.28);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{0.5}}, {1}, {1})[0], 0.96);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForFirstDerivativeAtLastKnot) { // NOLINT
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{1.0}}, {0}, {1})[0], -1.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{1.0}}, {1}, {1})[0], 0.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForSecondDerivativeAtFirstKnot) { // NOLINT
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{0.0}}, {0}, {2})[0], -4.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{0.0}}, {1}, {2})[0], 0.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForSecondDerivativeAtValueBetweenKnots) { // NOLINT
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{0.5}}, {0}, {2})[0],
              DoubleNear(-0.512, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{0.5}}, {1}, {2})[0],
              DoubleNear(-2.816, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForSecondDerivativeAtLastKnot) { // NOLINT
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{1.0}}, {0}, {2})[0], 1.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{1.0}}, {1}, {2})[0], -1.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForThirdDerivativeAtFirstKnot) { // NOLINT
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{0.0}}, {0}, {3})[0], 0.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{0.0}}, {1}, {3})[0], -12.0);
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForThirdDerivativeAtValueBetweenKnots) { // NOLINT
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{0.5}}, {0}, {3})[0],
              DoubleNear(7.3728, util::NumericSettings<double>::kEpsilon()));
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{0.5}}, {1}, {3})[0],
              DoubleNear(2.1504, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(NurbsDerivativeEx4_2, ReturnsCorrectValuesForThirdDerivativeAtLastKnot) { // NOLINT
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{1.0}}, {0}, {3})[0], 0.0);
  ASSERT_THAT(nurbs->EvaluateDerivative({ParamCoord{1.0}}, {1}, {3})[0], 3.0);
}

class ANURBSWithSplineGenerator : public Test {
 public:
  ANURBSWithSplineGenerator() {
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector =
        {std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                            ParamCoord{2}, ParamCoord{3},
                                                            ParamCoord{3}, ParamCoord{3}}))};
    std::array<Degree, 1> degree = {Degree{2}};
    std::vector<double> weights = {1, 4, 1, 1, 1};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({5.0, -1.0}))
    };
    spl::NURBSGenerator<1> nurbs_generator(knot_vector, degree, control_points, weights);
    nurbs = std::make_unique<spl::NURBS<1>>(nurbs_generator);
  }

 protected:
  std::unique_ptr<spl::NURBS<1>> nurbs;
};

TEST_F(ANURBSWithSplineGenerator, Returns1_4For1AndDim0) { // NOLINT
  ASSERT_THAT(nurbs->Evaluate({ParamCoord{1.0}}, {0})[0], DoubleNear(1.4, util::NumericSettings<double>::kEpsilon()));
}
