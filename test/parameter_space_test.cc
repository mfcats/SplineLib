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

#include "parameter_space.h"

using testing::Test;
using testing::DoubleEq;

class A1DParameterSpace : public Test {
 public:
  A1DParameterSpace() {
    std::array<baf::KnotVector, 1> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2}, ParamCoord{3},
                          ParamCoord{4}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}})};
    degree_ = {2};
    knot_vector_[0] = std::make_shared<baf::KnotVector>(knot_vector[0]);
    parameter_space = spl::ParameterSpace<1>(knot_vector_, degree_);
  }

 protected:
  spl::ParameterSpace<1> parameter_space;
  std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_;
  std::array<int, 1> degree_;
};

TEST_F(A1DParameterSpace, returnsCorrectDegree) {
  ASSERT_THAT(parameter_space.GetDegree(0), 2);
}

TEST_F(A1DParameterSpace, returns3_0ForFifthKnot) {
  ASSERT_THAT((*parameter_space.GetKnotVector(0))[5].get(), DoubleEq(3.0));
}

TEST_F(A1DParameterSpace, returnsCorrectBasisFunctionValuesForParamCoord0_5) {
  std::vector<double> values = {0.25, 0.625, 0.125, 0, 0, 0, 0, 0};
  for (int i = 0; i < 8; ++i) {
    ASSERT_THAT(parameter_space.GetBasisFunctions({i}, {ParamCoord(0.5)}), DoubleEq(values[i]));
  }
}

TEST_F(A1DParameterSpace, returnsCorrectBasisFunctionDerivativeValuesForParamCoord0_5AndDerivative1) {
  std::vector<double> values = {-1, 0.5, 0.5, 0, 0, 0, 0, 0};
  for (int i = 0; i < 8; ++i) {
    ASSERT_THAT(parameter_space.GetBasisFunctionDerivatives({i}, {ParamCoord(0.5)}, {1}), DoubleEq(values[i]));
  }
}

class A2DParameterSpace : public Test {
 public:
  A2DParameterSpace() {
    std::array<baf::KnotVector, 2> knot_vector =
        {baf::KnotVector({std::vector<ParamCoord>({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                   ParamCoord{1}, ParamCoord{1}})}),
         baf::KnotVector({std::vector<ParamCoord>({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2},
                                                   ParamCoord{3}, ParamCoord{3}})})};
    degree_ = {2, 1};
    knot_vector_[0] = std::make_shared<baf::KnotVector>(knot_vector[0]);
    knot_vector_[1] = std::make_shared<baf::KnotVector>(knot_vector[1]);
    parameter_space = spl::ParameterSpace<2>(knot_vector_, degree_);
  }

 protected:
  spl::ParameterSpace<2> parameter_space;
  std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector_;
  std::array<int, 2> degree_;
};

TEST_F(A2DParameterSpace, returnsDegree2ForFirstDimension) {
  ASSERT_THAT(parameter_space.GetDegree(0), 2);
}

TEST_F(A2DParameterSpace, returnsDegree1ForSecondDimension) {
  ASSERT_THAT(parameter_space.GetDegree(1), 1);
}

TEST_F(A2DParameterSpace, returns1_0ForFourthKnotOfFirstKnotVector) {
  ASSERT_THAT((*parameter_space.GetKnotVector(0))[3].get(), DoubleEq(1.0));
}

TEST_F(A2DParameterSpace, returns2_0ForFourthKnotOfSecondKnotVector) {
  ASSERT_THAT((*parameter_space.GetKnotVector(1))[3].get(), DoubleEq(2.0));
}

TEST_F(A2DParameterSpace, returnsCorrectBasisFunctionValue) {
  ASSERT_THAT(parameter_space.GetBasisFunctions({1, 1}, {ParamCoord(0.5), ParamCoord(0.25)}), DoubleEq(0.125));
}

TEST_F(A2DParameterSpace, returnsCorrectBasisFunctionDerivativeValue) {
  ASSERT_THAT(parameter_space.GetBasisFunctionDerivatives({1, 1}, {ParamCoord(0.5), ParamCoord(0.25)}, {0, 1}),
              DoubleEq(0.5));
}

class A3DParameterSpace : public Test {
 public:
  A3DParameterSpace() {
    std::array<baf::KnotVector, 3> knot_vector =
        {baf::KnotVector({std::vector<ParamCoord>({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                   ParamCoord{1}, ParamCoord{1}})}),
         baf::KnotVector({std::vector<ParamCoord>({ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{0.9}})}),
         baf::KnotVector({std::vector<ParamCoord>({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2},
                                                   ParamCoord{3}, ParamCoord{3}})})};
    degree_ = {2, 0, 1};
    knot_vector_[0] = std::make_shared<baf::KnotVector>(knot_vector[0]);
    knot_vector_[1] = std::make_shared<baf::KnotVector>(knot_vector[1]);
    knot_vector_[2] = std::make_shared<baf::KnotVector>(knot_vector[2]);
    parameter_space = spl::ParameterSpace<3>(knot_vector_, degree_);
  }

 protected:
  spl::ParameterSpace<3> parameter_space;
  std::array<std::shared_ptr<baf::KnotVector>, 3> knot_vector_;
  std::array<int, 3> degree_;
};

TEST_F(A3DParameterSpace, returnsDegree2ForFirstDimension) {
  ASSERT_THAT(parameter_space.GetDegree(0), 2);
}

TEST_F(A3DParameterSpace, returnsDegree0ForSecondDimension) {
  ASSERT_THAT(parameter_space.GetDegree(1), 0);
}

TEST_F(A3DParameterSpace, returnsDegree1ForThirdDimension) {
  ASSERT_THAT(parameter_space.GetDegree(2), 1);
}

TEST_F(A3DParameterSpace, returns1_0ForFourthKnotOfFirstKnotVector) {
  ASSERT_THAT((*parameter_space.GetKnotVector(0))[3].get(), DoubleEq(1.0));
}

TEST_F(A3DParameterSpace, returns0_3ForSecondKnotOfSecondKnotVector) {
  ASSERT_THAT((*parameter_space.GetKnotVector(1))[1].get(), DoubleEq(0.3));
}

TEST_F(A3DParameterSpace, returns2_0ForFourthKnotOfThirdKnotVector) {
  ASSERT_THAT((*parameter_space.GetKnotVector(2))[3].get(), DoubleEq(2.0));
}

TEST_F(A3DParameterSpace, returnsCorrectBasisFunctionValue) {
  ASSERT_THAT(parameter_space.GetBasisFunctions({1, 1, 1}, {ParamCoord(0.5), ParamCoord(0.5), ParamCoord(0.25)}),
              DoubleEq(0.125));
}

TEST_F(A3DParameterSpace, returnsCorrectBasisFunctionDerivativeValue) {
  ASSERT_THAT(parameter_space.GetBasisFunctionDerivatives({1, 1, 1},
                                                          {ParamCoord(0.5), ParamCoord(0.5), ParamCoord(0.25)},
                                                          {0, 0, 1}), DoubleEq(0.5));
}

