/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/util/numeric_operations.h"

using testing::Test;

using namespace splinelib::src;

class AnIntegerForNumericOperations : public Test {
 public:
  AnIntegerForNumericOperations() = default;

 protected:
  int an_integer_{10};
};

TEST_F(AnIntegerForNumericOperations, IsIncrementedTo11) {  // NOLINT
  an_integer_ = util::numeric_operations::increment(an_integer_);
  ASSERT_THAT(an_integer_, 11);
}

TEST_F(AnIntegerForNumericOperations, DecrementedTo9) {  // NOLINT
  an_integer_ = util::numeric_operations::decrement(an_integer_);
  ASSERT_THAT(an_integer_, 9);
}
