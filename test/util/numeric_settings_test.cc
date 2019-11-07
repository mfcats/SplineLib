/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/util/numeric_settings.h"

using testing::DoubleEq;
using testing::Eq;
using testing::Test;

using namespace splinelib::src;

class IntegerNumericSettings : public Test {
 public:
  IntegerNumericSettings() = default;
};

TEST_F(IntegerNumericSettings, ReturnCorrectIntegerEpsilon) {  // NOLINT
  ASSERT_THAT(util::numeric_settings::GetEpsilon<int>(), Eq(0));
}

TEST_F(IntegerNumericSettings, CompareTwoIntegersWithoutExplicitTolerance) {  // NOLINT
  ASSERT_THAT(util::numeric_settings::AreEqual<int>(0, 0), true);
  ASSERT_THAT(util::numeric_settings::AreEqual<int>(1, 0), false);
}

TEST_F(IntegerNumericSettings, CompareTwoIntegersWithTolerance1) {  // NOLINT
  ASSERT_THAT(util::numeric_settings::AreEqual<int>(9, 10, 1), true);
  ASSERT_THAT(util::numeric_settings::AreEqual<int>(8, 10, 1), false);
}

class DoubleNumericSettings : public Test {
 public:
  DoubleNumericSettings() = default;
};

TEST_F(DoubleNumericSettings, ReturnCorrectDoubleEpsilon) {  // NOLINT
  ASSERT_THAT(util::numeric_settings::GetEpsilon<double>(), DoubleEq(10 * std::numeric_limits<double>::epsilon()));
}

TEST_F(DoubleNumericSettings, CompareTwoDoublesWithoutExplicitTolerance) {  // NOLINT
  ASSERT_THAT(util::numeric_settings::AreEqual<double>(std::numeric_limits<double>::epsilon(), 0), true);
  ASSERT_THAT(util::numeric_settings::AreEqual<double>(11 * std::numeric_limits<double>::epsilon(), 0), false);
}

TEST_F(DoubleNumericSettings, CompareTwoDoublesWithTolerance0_00001) {  // NOLINT
  ASSERT_THAT(util::numeric_settings::AreEqual<double>(9.99999, 10, 1e-5), true);
  ASSERT_THAT(util::numeric_settings::AreEqual<double>(9.999989, 10, 1e-5), false);
}
