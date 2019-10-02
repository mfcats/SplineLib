/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/
#include <numeric>

#include "gmock/gmock.h"

#include "parameter_space.h"
#include "numeric_settings.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

using namespace splinelib::src;

class A1DParameterSpace : public Test {
 public:
  A1DParameterSpace() :
      degree_{Degree{2}},
      knot_vector_{std::make_shared<baf::KnotVector>(baf::KnotVector(
          {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1},
           ParametricCoordinate{2},
           ParametricCoordinate{3}, ParametricCoordinate{4}, ParametricCoordinate{4}, ParametricCoordinate{5},
           ParametricCoordinate{5},
           ParametricCoordinate{5}}))},
      parameter_space(spl::ParameterSpace<1>(knot_vector_, degree_)) {}

 protected:
  std::array<Degree, 1> degree_;
  baf::KnotVectors<1> knot_vector_;
  spl::ParameterSpace<1> parameter_space;
};

TEST_F(A1DParameterSpace, ThrowsForTooShortKnotVector) {  // NOLINT
  knot_vector_ = {std::make_shared<baf::KnotVector>(baf::KnotVector(
      {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{5},
       ParametricCoordinate{5}}))};
  ASSERT_THROW(spl::ParameterSpace<1>(knot_vector_, degree_), std::runtime_error);
}

TEST_F(A1DParameterSpace, ThrowsForMissingMultiplicityOfFirstKnot) {  // NOLINT
  knot_vector_ = {std::make_shared<baf::KnotVector>(baf::KnotVector(
      {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1}, ParametricCoordinate{5},
       ParametricCoordinate{5},
       ParametricCoordinate{5}}))};
  ASSERT_THROW(spl::ParameterSpace<1>(knot_vector_, degree_), std::runtime_error);
}

TEST_F(A1DParameterSpace, ThrowsForMissingMultiplicityOfLastKnot) {  // NOLINT
  knot_vector_ = {std::make_shared<baf::KnotVector>(baf::KnotVector(
      {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{4},
       ParametricCoordinate{5},
       ParametricCoordinate{5}}))};
  ASSERT_THROW(spl::ParameterSpace<1>(knot_vector_, degree_), std::runtime_error);
}

TEST_F(A1DParameterSpace, ReturnsCorrectDegree) {  // NOLINT
  ASSERT_THAT(parameter_space.GetDegree(0), Degree{2});
}

TEST_F(A1DParameterSpace, Returns3_0ForFifthKnot) {  // NOLINT
  ASSERT_THAT((*parameter_space.GetKnotVector(0))[5].get(), DoubleEq(3.0));
}

TEST_F(A1DParameterSpace, ReturnsCorrectBasisFunctionValuesForParametricCoordinate0_5) {  // NOLINT
  std::vector<double> values = {0.25, 0.625, 0.125, 0, 0, 0, 0, 0};
  for (int i = 0; i < 8; ++i) {
    ASSERT_THAT(parameter_space.GetBasisFunctions({i}, {ParametricCoordinate(0.5)}), DoubleEq(values[i]));
  }
}

TEST_F(A1DParameterSpace,
       ReturnsCorrectBasisFunctionDerivativeValuesForParametricCoordinate0_5AndDerivative1) {  // NOLINT
  std::vector<double> values = {-1, 0.5, 0.5, 0, 0, 0, 0, 0};
  for (int i = 0; i < 8; ++i) {
    ASSERT_THAT(parameter_space.GetBasisFunctionDerivatives({i}, {ParametricCoordinate(0.5)}, {1}),
                DoubleEq(values[i]));
  }
}

TEST_F(A1DParameterSpace, Returns1AsFirstNonZeroBasisFunctionsForParametricCoordinate1_5) { // NOLINT
  ASSERT_THAT(parameter_space.GetArrayOfFirstNonZeroBasisFunctions({ParametricCoordinate(1.5)})[0], 1);
}

TEST_F(A1DParameterSpace, Returns5AsFirstNonZeroBasisFunctionsForParametricCoordinate4) { // NOLINT
  ASSERT_THAT(parameter_space.GetArrayOfFirstNonZeroBasisFunctions({ParametricCoordinate(4.0)})[0], 5);
}

TEST_F(A1DParameterSpace, ThrowsIfParametricCoordinateOutsideKnotVectorRange) { // NOLINT
  ASSERT_THROW(parameter_space.ThrowIfParametricCoordinateOutsideKnotVectorRange({ParametricCoordinate(5.1)}),
               std::range_error);
}

TEST_F(A1DParameterSpace, InsertsKnotCorrectly) {  // NOLINT
  parameter_space.InsertKnot(ParametricCoordinate(0.5), 0);
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetNumberOfKnots(), 12);
  std::vector<double> values = {0, 0.5, 0.5, 0, 0, 0, 0, 0, 0};
  for (int i = 0; i < 9; ++i) {
    ASSERT_THAT(parameter_space.GetBasisFunctions({i}, {ParametricCoordinate(0.5)}), DoubleEq(values[i]));
  }
}

TEST_F(A1DParameterSpace, RemovesKnotCorrectly) {  // NOLINT
  parameter_space.RemoveKnot(ParametricCoordinate(1), 0);
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetNumberOfKnots(), 10);
  std::vector<double> values = {0.5625, 0.395833, 0.0416667, 0, 0, 0, 0};
  for (int i = 0; i < 7; ++i) {
    ASSERT_THAT(parameter_space.GetBasisFunctions({i}, {ParametricCoordinate(0.5)}), DoubleNear(values[i], 1e-6));
  }
}

TEST_F(A1DParameterSpace, EvaluatesThatCopiedSpaceEqualsOriginalSpace) {  // NOLINT
  spl::ParameterSpace<1> copy(parameter_space);
  ASSERT_THAT(parameter_space.AreEqual(copy), true);
}

class A2DParameterSpace : public Test {
 public:
  A2DParameterSpace() : degree_{Degree{2}, Degree{1}},
                        knot_vector_{
                            std::make_shared<baf::KnotVector>(baf::KnotVector{
                                ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0},
                                ParametricCoordinate{1},
                                ParametricCoordinate{1}, ParametricCoordinate{1}}),
                            std::make_shared<baf::KnotVector>(baf::KnotVector{
                                ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1},
                                ParametricCoordinate{2},
                                ParametricCoordinate{3}, ParametricCoordinate{3}})},
                        parameter_space(knot_vector_, degree_) {}

 protected:
  std::array<Degree, 2> degree_;
  baf::KnotVectors<2> knot_vector_;
  spl::ParameterSpace<2> parameter_space;
};

TEST_F(A2DParameterSpace, ReturnsDegree2ForFirstDimension) {  // NOLINT
  ASSERT_THAT(parameter_space.GetDegree(0), Degree{2});
}

TEST_F(A2DParameterSpace, ReturnsDegree1ForSecondDimension) {  // NOLINT
  ASSERT_THAT(parameter_space.GetDegree(1), Degree{1});
}

TEST_F(A2DParameterSpace, Returns1_0ForFourthKnotOfFirstKnotVector) {  // NOLINT
  ASSERT_THAT((*parameter_space.GetKnotVector(0))[3].get(), DoubleEq(1.0));
}

TEST_F(A2DParameterSpace, Returns2_0ForFourthKnotOfSecondKnotVector) {  // NOLINT
  ASSERT_THAT((*parameter_space.GetKnotVector(1))[3].get(), DoubleEq(2.0));
}

TEST_F(A2DParameterSpace, ReturnsCorrectBasisFunctionValue) {  // NOLINT
  ASSERT_THAT(parameter_space.GetBasisFunctions({1, 1}, {ParametricCoordinate(0.5), ParametricCoordinate(0.25)}),
              DoubleEq(0.125));
}

TEST_F(A2DParameterSpace, ReturnsCorrectBasisFunctionDerivativeValue) {  // NOLINT
  ASSERT_THAT(parameter_space.GetBasisFunctionDerivatives({1, 1},
                                                          {ParametricCoordinate(0.5), ParametricCoordinate(0.25)},
                                                          {0, 1}), DoubleEq(0.5));
}

TEST_F(A2DParameterSpace, EvaluatesThatCopiedSpaceEqualsOriginalSpace) {  // NOLINT
  spl::ParameterSpace<2> copy(parameter_space);
  ASSERT_THAT(parameter_space.AreEqual(copy), true);
}

class A3DParameterSpace : public Test {
 public:
  A3DParameterSpace() : degree_{Degree{2}, Degree{0}, Degree{1}},
                        knot_vector_{
                            std::make_shared<baf::KnotVector>(baf::KnotVector(
                                {std::vector<ParametricCoordinate>(
                                    {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0},
                                     ParametricCoordinate{1},
                                     ParametricCoordinate{1}, ParametricCoordinate{1}})})),
                            std::make_shared<baf::KnotVector>(baf::KnotVector(
                                {std::vector<ParametricCoordinate>(
                                    {ParametricCoordinate{0}, ParametricCoordinate{0.3}, ParametricCoordinate{0.6},
                                     ParametricCoordinate{0.9}})})),
                            std::make_shared<baf::KnotVector>(baf::KnotVector(
                                {std::vector<ParametricCoordinate>(
                                    {ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{1},
                                     ParametricCoordinate{2},
                                     ParametricCoordinate{3}, ParametricCoordinate{3}})}))},
                        parameter_space(knot_vector_, degree_) {}

 protected:
  std::array<Degree, 3> degree_;
  baf::KnotVectors<3> knot_vector_;
  spl::ParameterSpace<3> parameter_space;
};

TEST_F(A3DParameterSpace, ReturnsDegree2ForFirstDimension) {  // NOLINT
  ASSERT_THAT(parameter_space.GetDegree(0), Degree{2});
}

TEST_F(A3DParameterSpace, ReturnsDegree0ForSecondDimension) {  // NOLINT
  ASSERT_THAT(parameter_space.GetDegree(1), Degree{0});
}

TEST_F(A3DParameterSpace, ReturnsDegree1ForThirdDimension) {  // NOLINT
  ASSERT_THAT(parameter_space.GetDegree(2), Degree{1});
}

TEST_F(A3DParameterSpace, Returns1_0ForFourthKnotOfFirstKnotVector) {  // NOLINT
  ASSERT_THAT((*parameter_space.GetKnotVector(0))[3].get(), DoubleEq(1.0));
}

TEST_F(A3DParameterSpace, Returns0_3ForSecondKnotOfSecondKnotVector) {  // NOLINT
  ASSERT_THAT((*parameter_space.GetKnotVector(1))[1].get(), DoubleEq(0.3));
}

TEST_F(A3DParameterSpace, Returns2_0ForFourthKnotOfThirdKnotVector) {  // NOLINT
  ASSERT_THAT((*parameter_space.GetKnotVector(2))[3].get(), DoubleEq(2.0));
}

TEST_F(A3DParameterSpace, Returns3_0ForKnotVectorRangeThirdDimension) { // NOLINT
  ASSERT_THAT(parameter_space.GetKnotVectorRange(2), DoubleEq(3.0));
}

TEST_F(A3DParameterSpace, ReturnsCorrectBasisFunctionValue) {  // NOLINT
  ASSERT_THAT(parameter_space.GetBasisFunctions({1, 1, 1},
                                                {ParametricCoordinate(0.5), ParametricCoordinate(0.5),
                                                 ParametricCoordinate(0.25)}), DoubleEq(0.125));
}

TEST_F(A3DParameterSpace, ReturnsCorrectBasisFunctionDerivativeValue) {  // NOLINT
  ASSERT_THAT(parameter_space.GetBasisFunctionDerivatives({1, 1, 1},
                                                          {ParametricCoordinate(0.5), ParametricCoordinate(0.5),
                                                           ParametricCoordinate(0.25)}, {0, 0, 1}), DoubleEq(0.5));
}

TEST_F(A3DParameterSpace, EvaluatesThatCopiedSpaceEqualsOriginalSpace) {  // NOLINT
  spl::ParameterSpace<3> copy(parameter_space);
  ASSERT_THAT(parameter_space.AreEqual(copy), true);
}

