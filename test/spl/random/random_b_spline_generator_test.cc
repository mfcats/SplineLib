/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/spl/b_spline.h"
#include "test/spl/random/random_b_spline_generator.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

using namespace splinelib::src;
using namespace splinelib::test::spl::random;

class A1DRandomBSplineGenerator : public Test {  // NOLINT
 public:
  A1DRandomBSplineGenerator() : b_spline_generator_(RandomBSplineGenerator<1>(limits_, max_degree_, dimension_)) {}

 protected:
  const std::array<ParametricCoordinate, 2> limits_{{ParametricCoordinate{0.5}, ParametricCoordinate{2.75}}};
  const int max_degree_{10};
  const int dimension_{2};
  RandomBSplineGenerator<1> b_spline_generator_;
};

TEST_F(A1DRandomBSplineGenerator, CreatesCorrectKnots) {  // NOLINT
  spl::ParameterSpace<1> parameter_space = *b_spline_generator_.GetParameterSpace();
  size_t degree = static_cast<size_t>(parameter_space.GetDegree(0).Get());
  baf::KnotVector knot_vector = *parameter_space.GetKnotVector(0);
  size_t number_of_knots = knot_vector.GetNumberOfKnots();
  size_t i = 0;
  for (; i < degree + 1; ++i) {
    ASSERT_THAT(knot_vector[i].Get(), DoubleEq(limits_[0].Get()));
  }
  for (; i < number_of_knots - degree - 1; ++i) {
    ASSERT_THAT(knot_vector[i].Get() - knot_vector[i - 1].Get(),
                DoubleNear((limits_[1].Get() - limits_[0].Get()) / (number_of_knots - 2 * degree - 1), 0.00000001));
  }
  for (; i < number_of_knots; ++i) {
    ASSERT_THAT(knot_vector[i].Get(), DoubleEq(limits_[1].Get()));
  }
}

TEST_F(A1DRandomBSplineGenerator, RegardsMaximalDegree) {  // NOLINT
  ASSERT_LE(b_spline_generator_.GetParameterSpace()->GetDegree(0).Get(), max_degree_);
}

TEST_F(A1DRandomBSplineGenerator, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(b_spline_generator_.GetPhysicalSpace()->GetDimension(), dimension_);
}

class A2DRandomBSplineGenerator : public Test {  // NOLINT
 public:
  A2DRandomBSplineGenerator() : b_spline_generator_(RandomBSplineGenerator<2>(limits_, max_degree_, dimension_)) {}

 protected:
  const std::array<ParametricCoordinate, 2> limits_{{ParametricCoordinate{0.0}, ParametricCoordinate{10.0}}};
  const int max_degree_{8};
  const int dimension_{3};
  RandomBSplineGenerator<2> b_spline_generator_;
};

TEST_F(A2DRandomBSplineGenerator, RegardsParametricCoordinateLimits) {  // NOLINT
  spl::ParameterSpace<2> parameter_space = *b_spline_generator_.GetParameterSpace();
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
}

TEST_F(A2DRandomBSplineGenerator, RegardsMaximalDegree) {  // NOLINT
  spl::ParameterSpace<2> parameter_space = *b_spline_generator_.GetParameterSpace();
  ASSERT_LE(parameter_space.GetDegree(0).Get(), max_degree_);
  ASSERT_LE(parameter_space.GetDegree(1).Get(), max_degree_);
}

TEST_F(A2DRandomBSplineGenerator, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(b_spline_generator_.GetPhysicalSpace()->GetDimension(), dimension_);
}

class A3DRandomBSplineGenerator : public Test {  // NOLINT
 public:
  A3DRandomBSplineGenerator() : b_spline_generator_(RandomBSplineGenerator<3>(limits_, max_degree_, dimension_)) {}

 protected:
  const std::array<ParametricCoordinate, 2> limits_{{ParametricCoordinate{2.0}, ParametricCoordinate{3.0}}};
  const int max_degree_{20};
  const int dimension_{4};
  RandomBSplineGenerator<3> b_spline_generator_;
};

TEST_F(A3DRandomBSplineGenerator, RegardsParametricCoordinateLimits) {  // NOLINT
  spl::ParameterSpace<3> parameter_space = *b_spline_generator_.GetParameterSpace();
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(0)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(1)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(2)->GetFirstKnot().Get(), DoubleEq(limits_[0].Get()));
  ASSERT_THAT(parameter_space.GetKnotVector(2)->GetLastKnot().Get(), DoubleEq(limits_[1].Get()));
}

TEST_F(A3DRandomBSplineGenerator, RegardsMaximalDegree) {  // NOLINT
  spl::ParameterSpace<3> parameter_space = *b_spline_generator_.GetParameterSpace();
  ASSERT_LE(parameter_space.GetDegree(0).Get(), max_degree_);
  ASSERT_LE(parameter_space.GetDegree(1).Get(), max_degree_);
  ASSERT_LE(parameter_space.GetDegree(2).Get(), max_degree_);
}

TEST_F(A3DRandomBSplineGenerator, RegardsPhysicalDimension) {  // NOLINT
  ASSERT_LE(b_spline_generator_.GetPhysicalSpace()->GetDimension(), dimension_);
}
