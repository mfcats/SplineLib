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
#include "random_b_spline_generator.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class A1DRandomBSplineGenerator : public Test {  // NOLINT
 public:
  A1DRandomBSplineGenerator() {
    b_spline_generator_ = spl::RandomBSplineGenerator<1>(limits_, max_degree_, dimension_);
  }

 protected:
  spl::RandomBSplineGenerator<1> b_spline_generator_;
  const std::array<ParamCoord, 2> limits_{{ParamCoord{0.5}, ParamCoord{2.75}}};
  const int max_degree_{10};
  const int dimension_{2};
};

TEST_F(A1DRandomBSplineGenerator, CreatesCorrectKnots) {  // NOLINT
  spl::ParameterSpace<1> parameter_space = *b_spline_generator_.GetParameterSpace();
  size_t degree = static_cast<size_t>(parameter_space.GetDegree(0).get());
  baf::KnotVector knot_vector = *parameter_space.GetKnotVector(0);
  size_t i = 0;
  for (; i < degree + 1; ++i) {
    ASSERT_THAT(knot_vector.GetKnot(i).get(), DoubleEq(limits_[0].get()));
  }
  for (; i < knot_vector.GetNumberOfKnots() - degree - 1; ++i) {
    ASSERT_THAT(knot_vector.GetKnot(i).get() - knot_vector.GetKnot(i - 1).get(),
                DoubleNear((limits_[1].get() - limits_[0].get()) / (knot_vector.GetNumberOfKnots() - 2 * degree - 1),
                           0.00000001));
  }
  for (; i < knot_vector.GetNumberOfKnots(); ++i) {
    ASSERT_THAT(knot_vector.GetKnot(i).get(), DoubleEq(limits_[1].get()));
  }
}

TEST_F(A1DRandomBSplineGenerator, RegardsMaximalDegree) {  // NOLINT
  ASSERT_LE(b_spline_generator_.GetParameterSpace()->GetDegree(0).get(), max_degree_);
}

TEST_F(A1DRandomBSplineGenerator, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(b_spline_generator_.GetPhysicalSpace()->GetDimension(), dimension_);
}

class A2DRandomBSplineGenerator : public Test {  // NOLINT
 public:
  A2DRandomBSplineGenerator() {
    b_spline_generator_ = spl::RandomBSplineGenerator<2>(limits_, max_degree_, dimension_);
  }

 protected:
  spl::RandomBSplineGenerator<2> b_spline_generator_;
  const std::array<ParamCoord, 2> limits_{{ParamCoord{0.0}, ParamCoord{10.0}}};
  const int max_degree_{8};
  const int dimension_{3};
};

TEST_F(A2DRandomBSplineGenerator, RegardsParametricCoordinateLimits) {  // NOLINT
  spl::ParameterSpace<2> parameter_space = *b_spline_generator_.GetParameterSpace();
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetKnot(0).get(), DoubleEq(limits_[0].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetLastKnot().get(), DoubleEq(limits_[1].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetKnot(0).get(), DoubleEq(limits_[0].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetLastKnot().get(), DoubleEq(limits_[1].get()));
}

TEST_F(A2DRandomBSplineGenerator, RegardsMaximalDegree) {  // NOLINT
  spl::ParameterSpace<2> parameter_space = *b_spline_generator_.GetParameterSpace();
  ASSERT_LE(parameter_space.GetDegree(0).get(), max_degree_);
  ASSERT_LE(parameter_space.GetDegree(1).get(), max_degree_);
}

TEST_F(A2DRandomBSplineGenerator, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(b_spline_generator_.GetPhysicalSpace()->GetDimension(), dimension_);
}

class A3DRandomBSplineGenerator : public Test {  // NOLINT
 public:
  A3DRandomBSplineGenerator() {
    b_spline_generator_ = spl::RandomBSplineGenerator<3>(limits_, max_degree_, dimension_);
  }

 protected:
  spl::RandomBSplineGenerator<3> b_spline_generator_;
  const std::array<ParamCoord, 2> limits_{{ParamCoord{2.0}, ParamCoord{3.0}}};
  const int max_degree_{20};
  const int dimension_{4};
};

TEST_F(A3DRandomBSplineGenerator, RegardsParametricCoordinateLimits) {  // NOLINT
  spl::ParameterSpace<3> parameter_space = *b_spline_generator_.GetParameterSpace();
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetKnot(0).get(), DoubleEq(limits_[0].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetLastKnot().get(), DoubleEq(limits_[1].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetKnot(0).get(), DoubleEq(limits_[0].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetLastKnot().get(), DoubleEq(limits_[1].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(2)->GetKnot(0).get(), DoubleEq(limits_[0].get()));
  ASSERT_THAT(parameter_space.GetKnotVector(2)->GetLastKnot().get(), DoubleEq(limits_[1].get()));
}

TEST_F(A3DRandomBSplineGenerator, RegardsMaximalDegree) {  // NOLINT
  spl::ParameterSpace<3> parameter_space = *b_spline_generator_.GetParameterSpace();
  ASSERT_LE(parameter_space.GetDegree(0).get(), max_degree_);
  ASSERT_LE(parameter_space.GetDegree(1).get(), max_degree_);
  ASSERT_LE(parameter_space.GetDegree(2).get(), max_degree_);
}

TEST_F(A3DRandomBSplineGenerator, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(b_spline_generator_.GetPhysicalSpace()->GetDimension(), dimension_);
}
