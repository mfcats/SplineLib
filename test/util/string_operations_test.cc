/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/util/string_operations.h"

using testing::Test;
using testing::DoubleEq;

using namespace splinelib::src;

class StringOperations : public Test {
 public:
  StringOperations() = default;
};

TEST_F(StringOperations, CheckIfStringStartsWithhello) {  // NOLINT
  ASSERT_THAT(util::string_operations::StartsWith("hello world", "hello"), true);
  ASSERT_THAT(util::string_operations::StartsWith("Hello World", "hello"), false);
  ASSERT_THAT(util::string_operations::StartsWith("Oh hello", "hello"), false);
  ASSERT_THAT(util::string_operations::StartsWith("hella", "hello"), false);
  ASSERT_THAT(util::string_operations::StartsWith("hell", "hello"), false);
}

TEST_F(StringOperations, CheckIfStringEndsWithhello) {  // NOLINT
  ASSERT_THAT(util::string_operations::EndsWith("hello world", "hello"), false);
  ASSERT_THAT(util::string_operations::EndsWith("Oh Hello", "hello"), false);
  ASSERT_THAT(util::string_operations::EndsWith("Oh hello", "hello"), true);
  ASSERT_THAT(util::string_operations::EndsWith("bello", "hello"), false);
  ASSERT_THAT(util::string_operations::EndsWith("ello", "hello"), false);
}

TEST_F(StringOperations, SplitStringDelimitedBySpaces) {  // NOLINT
  std::string undivided = "a string to be splitted at spaces";
  std::vector<std::string> divided = util::string_operations::SplitStringAtDelimiter(undivided, ' ');
  ASSERT_THAT(divided.size(), 7);
  ASSERT_THAT(divided[3], "be");
}

TEST_F(StringOperations, TrimStringFromSpacesAndSquareBrackets) {  // NOLINT
  std::string untrimmed = "[ 7.88  ]   ";
  ASSERT_THAT(util::string_operations::TrimSpacesAndSquareBrackets(untrimmed), "7.88");
}

TEST_F(StringOperations, ConvertStringToDoubleForConvertibleInputFormat) {  // NOLINT
  ASSERT_THROW(util::string_operations::ConvertStringToNumber<double>("8,97"), std::invalid_argument);
  ASSERT_THAT(util::string_operations::ConvertStringToNumber<double>("8.97"), DoubleEq(8.97));
  ASSERT_THAT(util::string_operations::ConvertStringToNumber<double>("8"), DoubleEq(8));
  ASSERT_THAT(util::string_operations::ConvertStringToNumber<double>("8.0e-5"), DoubleEq(0.00008));
  ASSERT_THAT(util::string_operations::ConvertStringToNumber<double>("8e-5"), DoubleEq(0.00008));
  ASSERT_THAT(util::string_operations::ConvertStringToNumber<double>("8.0E-5"), DoubleEq(0.00008));
  ASSERT_THAT(util::string_operations::ConvertStringToNumber<double>("8E-5"), DoubleEq(0.00008));
}

TEST_F(StringOperations, ConvertStringToIntForConvertibleInputFormat) {  // NOLINT
  ASSERT_THROW(util::string_operations::ConvertStringToNumber<int>("8,97"), std::invalid_argument);
  ASSERT_THAT(util::string_operations::ConvertStringToNumber<int>("8.97"), 8);
  ASSERT_THAT(util::string_operations::ConvertStringToNumber<int>("8"), 8);
  ASSERT_THAT(util::string_operations::ConvertStringToNumber<int>("8.0e-5"), 0);
}

TEST_F(StringOperations, ConvertStringToBoolForConvertibleInputFormat) {  // NOLINT
  ASSERT_THROW(util::string_operations::ConvertStringToNumber<bool>("8,97"), std::invalid_argument);
  ASSERT_THAT(util::string_operations::ConvertStringToNumber<bool>("8.97"), true);
  ASSERT_THAT(util::string_operations::ConvertStringToNumber<bool>("0"), false);
}

TEST_F(StringOperations, ConvertStringVectorToIntVector) {  // NOLINT
  std::vector<std::string> string_vector = {"6.7", "9.9", "7.77", "1.1e-3"};
  std::vector<int> int_vector = util::string_operations::ConvertStringVectorToNumberVector<int>(string_vector);
  ASSERT_THAT(string_vector.size(), int_vector.size());
  ASSERT_THAT(int_vector.front(), 6);
  ASSERT_THAT(int_vector.back(), 0);
}

TEST_F(StringOperations, ConvertStringVectorToDoubleVector) {  // NOLINT
  std::vector<std::string> string_vector = {"6.7", "9.9", "7.77", "1.1e-3"};
  std::vector<double> double_vector = util::string_operations::ConvertStringVectorToNumberVector<double>(string_vector);
  ASSERT_THAT(string_vector.size(), double_vector.size());
  ASSERT_THAT(double_vector.front(), DoubleEq(6.7));
  ASSERT_THAT(double_vector.back(), DoubleEq(0.0011));
}

TEST_F(StringOperations, SplitStringDelimitedByCommaOrSemicolonIntoIntVector) {  // NOLINT
  std::string undivided = "6.0,5.5;7.0,7,10.10;0.99,1,8.8,9;8";
  std::vector<int> divided = util::string_operations::ConvertDelimitedStringToNumberVector<int>(undivided);
  ASSERT_THAT(divided.size(), 9);
  ASSERT_THAT(divided[3], 7);
  ASSERT_THAT(divided[4], 10);
  ASSERT_THAT(divided.back(), 9);
}

TEST_F(StringOperations, SplitStringDelimitedByCommaOrSemicolonIntoDoubleVector) {  // NOLINT
  std::string undivided = "6.0,5.5;7.0,7,10.10;0.99,1,8.8,9;8";
  std::vector<double> divided = util::string_operations::ConvertDelimitedStringToNumberVector<double>(undivided);
  ASSERT_THAT(divided.size(), 9);
  ASSERT_THAT(divided[3], DoubleEq(7.0));
  ASSERT_THAT(divided[4], DoubleEq(10.1));
  ASSERT_THAT(divided.back(), DoubleEq(9.0));
}
