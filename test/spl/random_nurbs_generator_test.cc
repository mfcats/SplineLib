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

#include "nurbs.h"
#include "random_nurbs_generator.h"

using testing::Test;
using testing::DoubleEq;

class A1DRandomNURBSGenerator : public Test {  // NOLINT
 public:
  A1DRandomNURBSGenerator() {
    nurbs_generator_ = spl::RandomNURBSGenerator<1>(limits_, max_degree_, dimension_);
  }

 protected:
  spl::RandomNURBSGenerator<1> nurbs_generator_;
  const std::array<ParamCoord, 2> limits_{{ParamCoord{0.5}, ParamCoord{2.75}}};
  const int max_degree_{10};
  const int dimension_{2};
};

TEST_F(A1DRandomNURBSGenerator, RegardsParametricCoordinateLimits) {  // NOLINT
  ASSERT_THAT(nurbs_generator_.GetParameterSpace()->GetKnotVector(0)->GetKnot(0).get(), DoubleEq(limits_[0].get()));
  ASSERT_THAT(nurbs_generator_.GetParameterSpace()->GetKnotVector(0)->GetLastKnot().get(),
              DoubleEq(limits_[1].get()));
}

TEST_F(A1DRandomNURBSGenerator, RegardsMaximalDegree) {  // NOLINT
  ASSERT_LE(nurbs_generator_.GetParameterSpace()->GetDegree(0).get(), max_degree_);
}

TEST_F(A1DRandomNURBSGenerator, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(nurbs_generator_.GetWeightedPhysicalSpace()->GetDimension(), dimension_);
}

class A2DRandomNURBSGenerator : public Test {  // NOLINT
 public:
  A2DRandomNURBSGenerator() {
    nurbs_generator_ = spl::RandomNURBSGenerator<2>(limits_, max_degree_, dimension_);
  }

 protected:
  spl::RandomNURBSGenerator<2> nurbs_generator_;
  const std::array<ParamCoord, 2> limits_{{ParamCoord{0.0}, ParamCoord{10.0}}};
  const int max_degree_{8};
  const int dimension_{3};
};

TEST_F(A2DRandomNURBSGenerator, RegardsParametricCoordinateLimits) {  // NOLINT
  spl::ParameterSpace<2> parameter_space = *nurbs_generator_.GetParameterSpace();
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetKnot(0).get(), DoubleEq(limits_[0].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetLastKnot().get(), DoubleEq(limits_[1].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetKnot(0).get(), DoubleEq(limits_[0].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetLastKnot().get(), DoubleEq(limits_[1].get()));
}

TEST_F(A2DRandomNURBSGenerator, RegardsMaximalDegree) {  // NOLINT
  spl::ParameterSpace<2> parameter_space = *nurbs_generator_.GetParameterSpace();
  ASSERT_LE(parameter_space.GetDegree(0).get(), max_degree_);
  ASSERT_LE(parameter_space.GetDegree(1).get(), max_degree_);
}

TEST_F(A2DRandomNURBSGenerator, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(nurbs_generator_.GetWeightedPhysicalSpace()->GetDimension(), dimension_);
}

class A3DRandomNURBSGenerator : public Test {  // NOLINT
 public:
  A3DRandomNURBSGenerator() {
    nurbs_generator_ = spl::RandomNURBSGenerator<3>(limits_, max_degree_, dimension_);
  }

 protected:
  spl::RandomNURBSGenerator<3> nurbs_generator_;
  const std::array<ParamCoord, 2> limits_{{ParamCoord{2.0}, ParamCoord{3.0}}};
  const int max_degree_{20};
  const int dimension_{4};
};

TEST_F(A3DRandomNURBSGenerator, RegardsParametricCoordinateLimits) {  // NOLINT
  spl::ParameterSpace<3> parameter_space = *nurbs_generator_.GetParameterSpace();
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetKnot(0).get(), DoubleEq(limits_[0].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetLastKnot().get(), DoubleEq(limits_[1].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetKnot(0).get(), DoubleEq(limits_[0].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetLastKnot().get(), DoubleEq(limits_[1].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(2)->GetKnot(0).get(), DoubleEq(limits_[0].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(2)->GetLastKnot().get(), DoubleEq(limits_[1].get()));
}

TEST_F(A3DRandomNURBSGenerator, RegardsMaximalDegree) {  // NOLINT
  spl::ParameterSpace<3> parameter_space = *nurbs_generator_.GetParameterSpace();
  ASSERT_LE(parameter_space.GetDegree(0).get(), max_degree_);
  ASSERT_LE(parameter_space.GetDegree(1).get(), max_degree_);
  ASSERT_LE(parameter_space.GetDegree(2).get(), max_degree_);
}

TEST_F(A3DRandomNURBSGenerator, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(nurbs_generator_.GetWeightedPhysicalSpace()->GetDimension(), dimension_);
}
