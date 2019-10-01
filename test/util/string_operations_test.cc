/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "string_operations.h"

#include "gmock/gmock.h"

using testing::Test;
using testing::DoubleEq;

using namespace splinelib::src;

class StringOperations : public Test {
 public:
  StringOperations() = default;
};

TEST_F(StringOperations, SplitStringDelimitedBySpaces) {  // NOLINT
  std::string undivided = "a string to be splitted at spaces";
  std::vector<std::string> divided = util::StringOperations::split(undivided, ' ');
  ASSERT_THAT(divided.size(), 7);
  ASSERT_THAT(divided[3], "be");
}

TEST_F(StringOperations, CheckIfStringStartsWithHello) {  // NOLINT
  ASSERT_THAT(util::StringOperations::StartsWith("hello world", "hello"), true);
  ASSERT_THAT(util::StringOperations::StartsWith("Oh hello", "hello"), false);
  ASSERT_THAT(util::StringOperations::StartsWith("hella", "hello"), false);
}

TEST_F(StringOperations, CheckIfStringEndsWithHello) {  // NOLINT
  ASSERT_THAT(util::StringOperations::EndsWith("hello world", "hello"), false);
  ASSERT_THAT(util::StringOperations::EndsWith("Oh hello", "hello"), true);
  ASSERT_THAT(util::StringOperations::EndsWith("bello", "hello"), false);
}

TEST_F(StringOperations, TrimString) {  // NOLINT
  std::string untrimmed = "[ 7.88  ] ";
  ASSERT_THAT(util::StringOperations::trim(untrimmed), "7.88");
}

TEST_F(StringOperations, ConvertStringToDouble) {  // NOLINT
  ASSERT_THROW(util::StringOperations::StringToDouble("8,97"), std::invalid_argument);
}

TEST_F(StringOperations, SplitStringDelimitedByCommaOrSemicolonInDoubleVector) {  // NOLINT
  std::string undivided = "6.0,5.5;7.0,7,10.10;0.99,1,8.8,9;";
  std::vector<double> divided = util::StringOperations::DelimitedStringToVector(undivided);
  ASSERT_THAT(divided.size(), 9);
  ASSERT_THAT(divided[3], DoubleEq(7.0));
  ASSERT_THAT(divided.back(), DoubleEq(9.0));
}

TEST_F(StringOperations, ConvertStringVectorToIntVector) {  // NOLINT
  std::vector<std::string> string_vector = {"6.7", "9.9", "7.77", "0.11"};
  std::vector<int> int_vector = util::StringOperations::StringVectorToNumberVector<int>(string_vector);
  ASSERT_THAT(string_vector.size(), int_vector.size());
  ASSERT_THAT(int_vector.front(), 6);
  ASSERT_THAT(int_vector.back(), 0);
}

TEST_F(StringOperations, ConvertStringVectorToDoubleVector) {  // NOLINT
  std::vector<std::string> string_vector = {"6.7", "9.9", "7.77", "0.11"};
  std::vector<double> double_vector = util::StringOperations::StringVectorToNumberVector<double>(string_vector);
  ASSERT_THAT(string_vector.size(), double_vector.size());
  ASSERT_THAT(double_vector.front(), DoubleEq(6.7));
  ASSERT_THAT(double_vector.back(), DoubleEq(0.11));
}
