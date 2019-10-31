/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include <array>
#include <vector>

#include "src/util/stl_container_access.h"

using testing::Test;
using testing::DoubleEq;

using namespace splinelib::src;

class AVectorForSTLContainerAccess : public Test {
 public:
  AVectorForSTLContainerAccess() : a_vector_({0.5, 1.2, 2.9, -3.4, -5}) {}

 protected:
  std::vector<double> a_vector_;
};

TEST_F(AVectorForSTLContainerAccess, Returns1_2ForConstSTLContainerAccessAtIndex1) {  // NOLINT
  ASSERT_THAT(GetValue(a_vector_, 1), DoubleEq(1.2));
}

TEST_F(AVectorForSTLContainerAccess, Returns1_2ForConstSTLContainerAccessWithNamedTypeIndex1) {  // NOLINT
  ASSERT_THAT(GetValue(a_vector_, Dimension{1}), DoubleEq(1.2));
}

TEST_F(AVectorForSTLContainerAccess, CanBeWrittenMinus10_1ToIndex1ViaSTLContainerAccess) {  // NOLINT
  GetValue(a_vector_, 1) = -10.1;
  ASSERT_THAT(GetValue(a_vector_, 1), DoubleEq(-10.1));
}

TEST_F(AVectorForSTLContainerAccess, CanBeWritten10_1ToNamedTypeIndex1ViaSTLContainerAccess) {  // NOLINT
  GetValue(a_vector_, Dimension{1}) = 10.1;
  ASSERT_THAT(GetValue(a_vector_, Dimension{1}), DoubleEq(10.1));
}

class AnArrayForSTLContainerAccess : public Test {
 public:
  AnArrayForSTLContainerAccess() : an_array_({0.5, 1.2, 2.9, -3.4, -5}) {}

 protected:
  std::array<double, 5> an_array_;
};

TEST_F(AnArrayForSTLContainerAccess, Returns1_2ForConstSTLContainerAccessAtIndex1) {  // NOLINT
  ASSERT_THAT(GetValue(an_array_, 1), DoubleEq(1.2));
}

TEST_F(AnArrayForSTLContainerAccess, Returns1_2ForConstSTLContainerAccessWithNamedTypeIndex1) {  // NOLINT
  ASSERT_THAT(GetValue(an_array_, Dimension{1}), DoubleEq(1.2));
}

TEST_F(AnArrayForSTLContainerAccess, CanBeWrittenMinus10_1ToIndex1ViaSTLContainerAccess) {  // NOLINT
  GetValue(an_array_, 1) = -10.1;
  ASSERT_THAT(GetValue(an_array_, 1), DoubleEq(-10.1));
}

TEST_F(AnArrayForSTLContainerAccess, CanBeWritten10_1ToNamedTypeIndex1ViaSTLContainerAccess) {  // NOLINT
  GetValue(an_array_, Dimension{1}) = 10.1;
  ASSERT_THAT(GetValue(an_array_, Dimension{1}), DoubleEq(10.1));
}
