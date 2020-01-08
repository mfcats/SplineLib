/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/spl/weighted_physical_space.h"

using testing::Test;
using testing::DoubleEq;

using namespace splinelib::src;

class A1DWeightedPhysicalSpace : public Test {
 public:
  A1DWeightedPhysicalSpace() {
    control_points = {
        spl::ControlPoint(std::vector<double>({0.0, 0.0})),
        spl::ControlPoint(std::vector<double>({1.0, 1.0})),
        spl::ControlPoint(std::vector<double>({3.0, 2.0})),
        spl::ControlPoint(std::vector<double>({4.0, 1.0})),
        spl::ControlPoint(std::vector<double>({5.0, -1.0}))
    };
    weights_ = {0.5, 0.75, 0.8, 1.0, 1.2};
    weighted_physical_space = spl::WeightedPhysicalSpace<1>(control_points, weights_, {5});
  }

 protected:
  spl::WeightedPhysicalSpace<1> weighted_physical_space;
  std::vector<spl::ControlPoint> control_points;
  std::vector<double> weights_;
};

TEST_F(A1DWeightedPhysicalSpace, ThrowsForDifferingNumberOfControlPointsAndWeights) {  // NOLINT
  control_points.emplace_back(std::vector<double>({0.0, 0.0}));
  ASSERT_THROW(spl::WeightedPhysicalSpace<1>(control_points, weights_, {6}), std::runtime_error);
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectWeight) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight(2).Get(), DoubleEq(0.8));
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectMinimumWeight) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetMinimumWeight(), DoubleEq(0.5));
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectFirstHomogenousControlPoint) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 1>{0})[Dimension{0}], DoubleEq(0.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 1>{0})[Dimension{1}], DoubleEq(0.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 1>{0})[Dimension{2}], DoubleEq(0.5));
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectInnerHomogenousControlPoint) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 1>{2})[Dimension{0}], DoubleEq(2.4));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 1>{2})[Dimension{1}], DoubleEq(1.6));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 1>{2})[Dimension{2}], DoubleEq(0.8));
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectLastHomogenousControlPoint) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 1>{4})[Dimension{0}], DoubleEq(6.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 1>{4})[Dimension{1}], DoubleEq(-1.2));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 1>{4})[Dimension{2}], DoubleEq(1.2));
}

TEST_F(A1DWeightedPhysicalSpace, AddsAndSetsNewControlPoint) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetTotalNumberOfControlPoints(), 5);
  // TODO(Corinna, Christoph): It should not be possible to add a control point without specifying a dimension.
  // weighted_physical_space.AddControlPoints(1);
  // weighted_physical_space.SetControlPoint(5, spl::ControlPoint(std::vector<double>({6.0, 0.0})));
  // weighted_physical_space.SetWeight(5, 4.3);
  // ASSERT_THAT(weighted_physical_space.GetTotalNumberOfControlPoints(), 6);
  // ASSERT_THAT(weighted_physical_space.GetControlPoint(5)[Dimension{0}], DoubleEq(6.0));
  // ASSERT_THAT(weighted_physical_space.GetWeight(5), DoubleEq(4.3));
  weighted_physical_space.SetControlPoint(4, spl::ControlPoint(std::vector<double>({6.0, 0.0})));
  weighted_physical_space.SetWeight(4, 4.3);
  ASSERT_THAT(weighted_physical_space.GetControlPoint(4)[Dimension{0}], DoubleEq(6.0));
  ASSERT_THAT(weighted_physical_space.GetWeight(4).Get(), DoubleEq(4.3));
}

TEST_F(A1DWeightedPhysicalSpace, RemovesControlPointAndWeight) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetTotalNumberOfControlPoints(), 5);
  weighted_physical_space.RemoveControlPoints(2);
  ASSERT_THAT(weighted_physical_space.GetTotalNumberOfControlPoints(), 3);
  ASSERT_THAT(weighted_physical_space.GetWeights().size(), 3);
}

TEST_F(A1DWeightedPhysicalSpace, ReturnsCorrectDividedControlPoints) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetDividedWeights(2, 2, 0).size(), 2);
  ASSERT_THAT(weighted_physical_space.GetDividedWeights(2, 2, 0)[0], DoubleEq(0.8));
}

TEST_F(A1DWeightedPhysicalSpace, EvaluatesThatCopiedSpaceEqualsOriginalSpace) {  // NOLINT
  spl::WeightedPhysicalSpace<1> copy(weighted_physical_space);
  ASSERT_THAT(weighted_physical_space.AreEqual(copy), true);
}

class A2DWeightedPhysicalSpace : public Test {
 public:
  A2DWeightedPhysicalSpace() {
    control_points = {
        spl::ControlPoint(std::vector<double>({0.0, 0.0})),
        spl::ControlPoint(std::vector<double>({1.0, 1.0})),
        spl::ControlPoint(std::vector<double>({3.0, 2.0})),
        spl::ControlPoint(std::vector<double>({0.0, 2.0})),
        spl::ControlPoint(std::vector<double>({1.5, 2.5})),
        spl::ControlPoint(std::vector<double>({5.0, 1.0}))
    };
    weights_ = {0.5, 0.75, 0.8, 1.0, 1.2, 3.8};
    weighted_physical_space = spl::WeightedPhysicalSpace<2>(control_points, weights_, {3, 2});
  }

 protected:
  spl::WeightedPhysicalSpace<2> weighted_physical_space;
  std::vector<spl::ControlPoint> control_points;
  std::vector<double> weights_;
};

TEST_F(A2DWeightedPhysicalSpace, ThrowsForDifferingNumberOfControlPointsAndWeights) {  // NOLINT
  weights_.emplace_back(1.0);
  ASSERT_THROW(spl::WeightedPhysicalSpace<2>(control_points, weights_, {3, 2}), std::runtime_error);
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectFirstWeightFor2DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight({0, 0}).Get(), DoubleEq(0.5));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectInnerWeightFor2DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight({0, 1}).Get(), DoubleEq(1.0));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectLastWeightFor2DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight({2, 1}).Get(), DoubleEq(3.8));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectFirstWeightFor1DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight(0).Get(), DoubleEq(0.5));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectInnerWeightFor1DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight(3).Get(), DoubleEq(1.0));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectLastWeightFor1DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetWeight(5).Get(), DoubleEq(3.8));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectMinimumWeight) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetMinimumWeight(), DoubleEq(0.5));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectFirstHomogenousControlPointFor2DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 2>{0, 0})[Dimension{0}], DoubleEq(0.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 2>{0, 0})[Dimension{1}], DoubleEq(0.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 2>{0, 0})[Dimension{2}], DoubleEq(0.5));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectInnerHomogenousControlPointFor2DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 2>{1, 1})[Dimension{0}], DoubleEq(1.8));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 2>{1, 1})[Dimension{1}], DoubleEq(3.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 2>{1, 1})[Dimension{2}], DoubleEq(1.2));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectLastHomogenousControlPointFor2DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 2>{2, 1})[Dimension{0}], DoubleEq(19.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 2>{2, 1})[Dimension{1}], DoubleEq(3.8));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(
      std::array<int, 2>{2, 1})[Dimension{2}], DoubleEq(3.8));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectFirstHomogenousControlPointFor1DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(0)[Dimension{0}], DoubleEq(0.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(0)[Dimension{1}], DoubleEq(0.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(0)[Dimension{2}], DoubleEq(0.5));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectInnerHomogenousControlPointFor1DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(4)[Dimension{0}], DoubleEq(1.8));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(4)[Dimension{1}], DoubleEq(3.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(4)[Dimension{2}], DoubleEq(1.2));
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectLastHomogenousControlPointFor1DIndex) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(5)[Dimension{0}], DoubleEq(19.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(5)[Dimension{1}], DoubleEq(3.8));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(5)[Dimension{2}], DoubleEq(3.8));
}

TEST_F(A2DWeightedPhysicalSpace, AddsAndSetsNewControlPoint) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetTotalNumberOfControlPoints(), 6);
  // TODO(Corinna, Christoph): It should not be possible to add a control point without specifying a dimension.
  // weighted_physical_space.AddControlPoints(1);
  // weighted_physical_space.SetControlPoint(6, spl::ControlPoint(std::vector<double>({6.0, 0.0})));
  // weighted_physical_space.SetWeight(6, 4.3);
  // ASSERT_THAT(weighted_physical_space.GetTotalNumberOfControlPoints(), 7);
  // ASSERT_THAT(weighted_physical_space.GetControlPoint(6)[Dimension{0}], DoubleEq(6.0));
  // ASSERT_THAT(weighted_physical_space.GetWeight(6), DoubleEq(4.3));
  weighted_physical_space.SetControlPoint(5, spl::ControlPoint(std::vector<double>({6.0, 0.0})));
  weighted_physical_space.SetWeight(5, 4.3);
  ASSERT_THAT(weighted_physical_space.GetControlPoint(5)[Dimension{0}], DoubleEq(6.0));
  ASSERT_THAT(weighted_physical_space.GetWeight(5).Get(), DoubleEq(4.3));
}

TEST_F(A2DWeightedPhysicalSpace, Removes2ControlPointsAndWeights) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetTotalNumberOfControlPoints(), 6);
  weighted_physical_space.RemoveControlPoints(2);
  ASSERT_THAT(weighted_physical_space.GetTotalNumberOfControlPoints(), 4);
  ASSERT_THAT(weighted_physical_space.GetWeights().size(), 4);
}

TEST_F(A2DWeightedPhysicalSpace, ReturnsCorrectDividedWeights) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetDividedWeights(1, 2, 0).size(), 4);
  ASSERT_THAT(weighted_physical_space.GetDividedWeights(1, 2, 0)[2], DoubleEq(1.2));
}

TEST_F(A2DWeightedPhysicalSpace, EvaluatesThatCopiedSpaceEqualsOriginalSpace) {  // NOLINT
  spl::WeightedPhysicalSpace<2> copy(weighted_physical_space);
  ASSERT_THAT(weighted_physical_space.AreEqual(copy), true);
}

/* 1-dimensional nurbs spline with following properties :
 * KnotVector = {0, 0, 0, 1, 2, 3, 3, 3}
 * ControlPoints = {{0, 0}, {1, 1}, {3, 2}, {4, 1}, {5, -1}}
 * Weights = {1, 4, 1, 1, 1}
*/

class A1DWeightedPhysicalSpace_a : public Test {
 public:
  A1DWeightedPhysicalSpace_a() {
    control_points = {
        spl::ControlPoint(std::vector<double>({0.0, 0.0})),
        spl::ControlPoint(std::vector<double>({1.0, 1.0})),
        spl::ControlPoint(std::vector<double>({3.0, 2.0})),
        spl::ControlPoint(std::vector<double>({4.0, 1.0})),
        spl::ControlPoint(std::vector<double>({5.0, -1.0}))
    };
    weights_ = {1, 4, 1, 1, 1};
    weighted_physical_space = spl::WeightedPhysicalSpace<1>(control_points, weights_, {5});
  }

 protected:
  spl::WeightedPhysicalSpace<1> weighted_physical_space;
  std::vector<spl::ControlPoint> control_points;
  std::vector<double> weights_;
};

TEST_F(A1DWeightedPhysicalSpace_a, ReturnsHomogeneousControlPoint) {  // NOLINT
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(1)[Dimension{0}], DoubleEq(4.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(1)[Dimension{1}], DoubleEq(4.0));
  ASSERT_THAT(weighted_physical_space.GetHomogenousControlPoint(1)[Dimension{2}], DoubleEq(4.0));
}
