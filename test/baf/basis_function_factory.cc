/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <memory>

#include "gmock/gmock.h"

#include "basis_function_factory.h"

using testing::Test;

using namespace splinelib::src;

class ABasisFunctionFactory : public Test {
 public:
  ABasisFunctionFactory() : degree_{Degree{-1}},
                            knot_vector_({ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0},
                                          ParametricCoordinate{1}, ParametricCoordinate{1},
                                          ParametricCoordinate{1}}),
                            start_of_support_{KnotSpan{4}} {}

 protected:
  Degree degree_;
  baf::KnotVector knot_vector_;
  KnotSpan start_of_support_;
  baf::BasisFunctionFactory basis_function_factory;
};

TEST_F(ABasisFunctionFactory, throwsError) { //NOLINT
  ASSERT_THROW(basis_function_factory.CreateDynamic(knot_vector_, start_of_support_, degree_), std::runtime_error);
}
