/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/spl/nurbs.h"
#include "test/spl/random/random_spline_utils.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

using namespace splinelib::src;
using namespace splinelib::test;

class A1DRandomNURBS : public Test {  // NOLINT
 public:
  A1DRandomNURBS() : random_nurbs_(RandomSplineUtils<1>::GenerateRandomNURBS(limits_, max_degree_, dimension_)) {}

 protected:
  const std::array<ParametricCoordinate, 2> limits_{{ParametricCoordinate{0.5}, ParametricCoordinate{2.75}}};
  const int max_degree_{10};
  const int dimension_{2};
  std::shared_ptr<spl::NURBS<1>> random_nurbs_;
};

TEST_F(A1DRandomNURBS, CreatesCorrectKnots) {  // NOLINT
  int degree = random_nurbs_->GetDegree(0).Get();
  auto knot_vector = random_nurbs_->GetKnotVector(0);
  int i = 0;
  for (; i < degree + 1; ++i) {
    ASSERT_THAT((*knot_vector)[i].Get(), DoubleEq(limits_[0].Get()));
  }
  for (; i < (*knot_vector).GetNumberOfKnots() - degree - 1; ++i) {
    ASSERT_THAT((*knot_vector)[i].Get() - (*knot_vector)[i - 1].Get(),
                DoubleNear((limits_[1].Get() - limits_[0].Get()) / ((*knot_vector).GetNumberOfKnots() - 2 * degree - 1),
                           0.00000001));
  }
  for (; i < (*knot_vector).GetNumberOfKnots(); ++i) {
    ASSERT_THAT((*knot_vector)[i].Get(), DoubleEq(limits_[1].Get()));
  }
}

TEST_F(A1DRandomNURBS, RegardsMaximalDegree) {  // NOLINT
  ASSERT_LE(random_nurbs_->GetDegree(0).Get(), max_degree_);
}

TEST_F(A1DRandomNURBS, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(random_nurbs_->GetPointDim(), dimension_);
}

class A2DRandomNURBS : public Test {  // NOLINT
 public:
  A2DRandomNURBS() : random_nurbs_(RandomSplineUtils<2>::GenerateRandomNURBS(limits_, max_degree_, dimension_)) {}

 protected:
  const std::array<ParametricCoordinate, 2> limits_{{ParametricCoordinate{0.0}, ParametricCoordinate{10.0}}};
  const int max_degree_{8};
  const int dimension_{3};
  std::shared_ptr<spl::NURBS<2>> random_nurbs_;
};

TEST_F(A2DRandomNURBS, RegardsParametricCoordinateLimits) {  // NOLINT
  ASSERT_THAT(random_nurbs_->GetKnotVector(0)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(random_nurbs_->GetKnotVector(0)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
  ASSERT_THAT(random_nurbs_->GetKnotVector(1)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(random_nurbs_->GetKnotVector(1)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
}

TEST_F(A2DRandomNURBS, RegardsMaximalDegree) {  // NOLINT
  ASSERT_LE(random_nurbs_->GetDegree(0).Get(), max_degree_);
  ASSERT_LE(random_nurbs_->GetDegree(1).Get(), max_degree_);
}

TEST_F(A2DRandomNURBS, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(random_nurbs_->GetPointDim(), dimension_);
}

class A3DRandomNURBS : public Test {  // NOLINT
 public:
  A3DRandomNURBS() : random_nurbs_(RandomSplineUtils<3>::GenerateRandomNURBS(limits_, max_degree_, dimension_)) {}

 protected:
  const std::array<ParametricCoordinate, 2> limits_{{ParametricCoordinate{2.0}, ParametricCoordinate{3.0}}};
  const int max_degree_{20};
  const int dimension_{4};
  std::shared_ptr<spl::NURBS<3>> random_nurbs_;
};

TEST_F(A3DRandomNURBS, RegardsParametricCoordinateLimits) {  // NOLINT
  ASSERT_THAT(random_nurbs_->GetKnotVector(0)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(random_nurbs_->GetKnotVector(0)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
  ASSERT_THAT(random_nurbs_->GetKnotVector(1)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(random_nurbs_->GetKnotVector(1)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
  ASSERT_THAT(random_nurbs_->GetKnotVector(2)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(random_nurbs_->GetKnotVector(2)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
}

TEST_F(A3DRandomNURBS, RegardsMaximalDegree) {  // NOLINT
  ASSERT_LE(random_nurbs_->GetDegree(0).Get(), max_degree_);
  ASSERT_LE(random_nurbs_->GetDegree(1).Get(), max_degree_);
  ASSERT_LE(random_nurbs_->GetDegree(2).Get(), max_degree_);
}

TEST_F(A3DRandomNURBS, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(random_nurbs_->GetPointDim(), dimension_);
}
