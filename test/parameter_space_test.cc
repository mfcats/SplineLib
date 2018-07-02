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
    knot_vector_ =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2}, ParamCoord{3},
                          ParamCoord{4}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}})};
    degree_ = {2};
    parameter_space = spl::ParameterSpace<1>(knot_vector_, degree_);
  }

 protected:
  spl::ParameterSpace<1> parameter_space;
  std::array<baf::KnotVector, 1> knot_vector_;
  std::array<int, 1> degree_;
};

TEST_F(A1DParameterSpace, returnsCorrectDegree) {
  ASSERT_THAT(parameter_space.GetDegree(0), 2);
}

