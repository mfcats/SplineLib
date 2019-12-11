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
#include "src/spl/random_nurbs_generator.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

using namespace splinelib::src;

class A1DRandomNURBSGenerator : public Test {  // NOLINT
 public:
  A1DRandomNURBSGenerator() : nurbs_generator_(spl::RandomNURBSGenerator<1>(limits_, max_degree_, dimension_)) {}

 protected:
  const std::array<ParametricCoordinate, 2> limits_{{ParametricCoordinate{0.5}, ParametricCoordinate{2.75}}};
  const int max_degree_{10};
  const int dimension_{2};
  spl::RandomNURBSGenerator<1> nurbs_generator_;
};

TEST_F(A1DRandomNURBSGenerator, CreatesCorrectKnots) {  // NOLINT
  spl::ParameterSpace<1> parameter_space = *nurbs_generator_.GetParameterSpace();
  int degree = parameter_space.GetDegree(0).Get();
  baf::KnotVector knot_vector = *parameter_space.GetKnotVector(0);
  int i = 0;
  for (; i < degree + 1; ++i) {
    ASSERT_THAT(knot_vector[i].Get(), DoubleEq(limits_[0].Get()));
  }
  for (; i < knot_vector.GetNumberOfKnots() - degree - 1; ++i) {
    ASSERT_THAT(knot_vector[i].Get() - knot_vector[i - 1].Get(),
                DoubleNear((limits_[1].Get() - limits_[0].Get()) / (knot_vector.GetNumberOfKnots() - 2 * degree - 1),
                           0.00000001));
  }
  for (; i < knot_vector.GetNumberOfKnots(); ++i) {
    ASSERT_THAT(knot_vector[i].Get(), DoubleEq(limits_[1].Get()));
  }
}

TEST_F(A1DRandomNURBSGenerator, RegardsMaximalDegree) {  // NOLINT
  ASSERT_LE(nurbs_generator_.GetParameterSpace()->GetDegree(0).Get(), max_degree_);
}

TEST_F(A1DRandomNURBSGenerator, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(nurbs_generator_.GetWeightedPhysicalSpace()->GetDimension(), dimension_);
}

class A2DRandomNURBSGenerator : public Test {  // NOLINT
 public:
  A2DRandomNURBSGenerator() : nurbs_generator_(spl::RandomNURBSGenerator<2>(limits_, max_degree_, dimension_)) {}

 protected:
  const std::array<ParametricCoordinate, 2> limits_{{ParametricCoordinate{0.0}, ParametricCoordinate{10.0}}};
  const int max_degree_{8};
  const int dimension_{3};
  spl::RandomNURBSGenerator<2> nurbs_generator_;
};

TEST_F(A2DRandomNURBSGenerator, RegardsParametricCoordinateLimits) {  // NOLINT
  spl::ParameterSpace<2> parameter_space = *nurbs_generator_.GetParameterSpace();
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
}

TEST_F(A2DRandomNURBSGenerator, RegardsMaximalDegree) {  // NOLINT
  spl::ParameterSpace<2> parameter_space = *nurbs_generator_.GetParameterSpace();
  ASSERT_LE(parameter_space.GetDegree(0).Get(), max_degree_);
  ASSERT_LE(parameter_space.GetDegree(1).Get(), max_degree_);
}

TEST_F(A2DRandomNURBSGenerator, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(nurbs_generator_.GetWeightedPhysicalSpace()->GetDimension(), dimension_);
}

class A3DRandomNURBSGenerator : public Test {  // NOLINT
 public:
  A3DRandomNURBSGenerator() : nurbs_generator_(spl::RandomNURBSGenerator<3>(limits_, max_degree_, dimension_)) {}

 protected:
  const std::array<ParametricCoordinate, 2> limits_{{ParametricCoordinate{2.0}, ParametricCoordinate{3.0}}};
  const int max_degree_{20};
  const int dimension_{4};
  spl::RandomNURBSGenerator<3> nurbs_generator_;
};

TEST_F(A3DRandomNURBSGenerator, RegardsParametricCoordinateLimits) {  // NOLINT
  spl::ParameterSpace<3> parameter_space = *nurbs_generator_.GetParameterSpace();
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(2)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(2)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
}

TEST_F(A3DRandomNURBSGenerator, RegardsMaximalDegree) {  // NOLINT
  spl::ParameterSpace<3> parameter_space = *nurbs_generator_.GetParameterSpace();
  ASSERT_LE(parameter_space.GetDegree(0).Get(), max_degree_);
  ASSERT_LE(parameter_space.GetDegree(1).Get(), max_degree_);
  ASSERT_LE(parameter_space.GetDegree(2).Get(), max_degree_);
}

TEST_F(A3DRandomNURBSGenerator, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(nurbs_generator_.GetWeightedPhysicalSpace()->GetDimension(), dimension_);
}
