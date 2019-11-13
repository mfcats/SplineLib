/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include <array>
#include <vector>

#include "gmock/gmock.h"

#include "src/util/stl_container_access.h"

using testing::DoubleEq;
using testing::Test;

using namespace splinelib::src;

// TODO(all): should there only be one ASSERT test per TEST_F instance or can some of the following tests be aggregated?

class AVectorForSTLContainerAccess : public Test {
 public:
  AVectorForSTLContainerAccess() : a_vector_({0.5, 1.2, 2.9, -3.4, -5}) {}

 protected:
  std::vector<double> a_vector_;
};

TEST_F(AVectorForSTLContainerAccess, Returns1_2ForConstSTLContainerAccessAtIndex1) {  // NOLINT
  ASSERT_THAT(GetValue(a_vector_, 1), DoubleEq(1.2));
}

TEST_F(AVectorForSTLContainerAccess, ThrowsExceptionForConstSTLContainerAccessAtIndex5OnlyInDebugMode) {  // NOLINT
#ifdef DEBUG
  ASSERT_THROW(GetValue(a_vector_, 5), std::out_of_range);
#else
  ASSERT_NO_THROW(GetValue(a_vector_, 5));
#endif
}

TEST_F(AVectorForSTLContainerAccess, CanBeWrittenMinus10_1ToIndex1ViaSTLContainerAccess) {  // NOLINT
  GetValue(a_vector_, 1) = -10.1;
  ASSERT_THAT(GetValue(a_vector_, 1), DoubleEq(-10.1));
}

TEST_F(AVectorForSTLContainerAccess, ThrowsExceptionForWritingToIndex5ViaSTLContainerAccessOnlyInDebugMode) {  // NOLINT
#ifdef DEBUG
  ASSERT_THROW(GetValue(a_vector_, 5) = 9.9, std::out_of_range);
#else
  ASSERT_NO_THROW(GetValue(a_vector_, 5) = 9.9);
#endif
}

TEST_F(AVectorForSTLContainerAccess, Returns1_2ForConstSTLContainerAccessWithNamedTypeIndex1) {  // NOLINT
  ASSERT_THAT(GetValue(a_vector_, Dimension{1}), DoubleEq(1.2));
}

TEST_F(AVectorForSTLContainerAccess,   // NOLINT
       ThrowsExceptionForConstSTLContainerAccessWithNamedTypeIndex5OnlyInDebugMode) {
#ifdef DEBUG
  ASSERT_THROW(GetValue(a_vector_, Dimension{5}), std::out_of_range);
#else
  ASSERT_NO_THROW(GetValue(a_vector_, Dimension{5}));
#endif
}

TEST_F(AVectorForSTLContainerAccess, CanBeWritten10_1ToNamedTypeIndex1ViaSTLContainerAccess) {  // NOLINT
  GetValue(a_vector_, Dimension{1}) = 10.1;
  ASSERT_THAT(GetValue(a_vector_, Dimension{1}), DoubleEq(10.1));
}

TEST_F(AVectorForSTLContainerAccess,  // NOLINT
       ThrowsExceptionForWritingToNamedTypeIndex5ViaSTLContainerAccessOnlyInDebugMode) {
#ifdef DEBUG
  ASSERT_THROW(GetValue(a_vector_, Dimension{5}) = 9.9, std::out_of_range);
#else
  ASSERT_NO_THROW(GetValue(a_vector_, Dimension{5}) = 9.9);
#endif
}

TEST_F(AVectorForSTLContainerAccess, Returns0_5ForConstSTLContainerAccessAtFront) {  // NOLINT
  ASSERT_THAT(GetFront(a_vector_), DoubleEq(0.5));
}

TEST_F(AVectorForSTLContainerAccess, ThrowsExceptionForEmptyConstSTLContainerAccessAtFrontOnlyInDebugMode) {  // NOLINT
  a_vector_.clear();
#ifdef DEBUG
  ASSERT_THROW(GetFront(a_vector_), std::out_of_range);
#else
  ASSERT_NO_THROW(GetFront(a_vector_));
#endif
}

TEST_F(AVectorForSTLContainerAccess, CanBeWritten10_1ToFrontViaSTLContainerAccess) {  // NOLINT
  GetFront(a_vector_) = 10.1;
  ASSERT_THAT(GetValue(a_vector_, 0), DoubleEq(10.1));
}

TEST_F(AVectorForSTLContainerAccess,  // NOLINT
       ThrowsExceptionForWritingToEmptyFrontViaSTLContainerAccessOnlyInDebugMode) {
  a_vector_.clear();
#ifdef DEBUG
  ASSERT_THROW(GetFront(a_vector_) = 9.9, std::out_of_range);
#else
  ASSERT_NO_THROW(GetFront(a_vector_) = 9.9);
#endif
}

TEST_F(AVectorForSTLContainerAccess, ReturnsMinus5ForConstSTLContainerAccessAtBack) {  // NOLINT
  ASSERT_THAT(GetBack(a_vector_), DoubleEq(-5));
}

TEST_F(AVectorForSTLContainerAccess, ThrowsExceptionForEmptyConstSTLContainerAccessAtBackOnlyInDebugMode) {  // NOLINT
  a_vector_.clear();
#ifdef DEBUG
  ASSERT_THROW(GetBack(a_vector_), std::out_of_range);
#else
  ASSERT_NO_THROW(GetBack(a_vector_));
#endif
}

TEST_F(AVectorForSTLContainerAccess, CanBeWritten10_1ToBackViaSTLContainerAccess) {  // NOLINT
  GetBack(a_vector_) = 10.1;
  ASSERT_THAT(GetValue(a_vector_, 4), DoubleEq(10.1));
}

TEST_F(AVectorForSTLContainerAccess,  // NOLINT
       ThrowsExceptionForWritingToEmptyBackViaSTLContainerAccessOnlyInDebugMode) {
  a_vector_.clear();
#ifdef DEBUG
  ASSERT_THROW(GetBack(a_vector_) = 9.9, std::out_of_range);
#else
  ASSERT_NO_THROW(GetBack(a_vector_) = 9.9);
#endif
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

TEST_F(AnArrayForSTLContainerAccess, ThrowsExceptionForConstSTLContainerAccessAtIndex5OnlyInDebugMode) {  // NOLINT
#ifdef DEBUG
  ASSERT_THROW(GetValue(an_array_, 5), std::out_of_range);
#else
  ASSERT_NO_THROW(GetValue(an_array_, 5));
#endif
}

TEST_F(AnArrayForSTLContainerAccess, CanBeWrittenMinus10_1ToIndex1ViaSTLContainerAccess) {  // NOLINT
  GetValue(an_array_, 1) = -10.1;
  ASSERT_THAT(GetValue(an_array_, 1), DoubleEq(-10.1));
}

TEST_F(AnArrayForSTLContainerAccess, ThrowsExceptionForWritingToIndex5ViaSTLContainerAccessOnlyInDebugMode) {  // NOLINT
#ifdef DEBUG
  ASSERT_THROW(GetValue(an_array_, 5) = 9.9, std::out_of_range);
#else
  ASSERT_NO_THROW(GetValue(an_array_, 5) = 9.9);
#endif
}

TEST_F(AnArrayForSTLContainerAccess, Returns1_2ForConstSTLContainerAccessWithNamedTypeIndex1) {  // NOLINT
  ASSERT_THAT(GetValue(an_array_, Dimension{1}), DoubleEq(1.2));
}

TEST_F(AnArrayForSTLContainerAccess,   // NOLINT
       ThrowsExceptionForConstSTLContainerAccessWithNamedTypeIndex5OnlyInDebugMode) {
#ifdef DEBUG
  ASSERT_THROW(GetValue(an_array_, Dimension{5}), std::out_of_range);
#else
  ASSERT_NO_THROW(GetValue(an_array_, Dimension{5}));
#endif
}

TEST_F(AnArrayForSTLContainerAccess, CanBeWritten10_1ToNamedTypeIndex1ViaSTLContainerAccess) {  // NOLINT
  GetValue(an_array_, Dimension{1}) = 10.1;
  ASSERT_THAT(GetValue(an_array_, Dimension{1}), DoubleEq(10.1));
}

TEST_F(AnArrayForSTLContainerAccess,  // NOLINT
       ThrowsExceptionForWritingToNamedTypeIndex5ViaSTLContainerAccessOnlyInDebugMode) {
#ifdef DEBUG
  ASSERT_THROW(GetValue(an_array_, Dimension{5}) = 9.9, std::out_of_range);
#else
  ASSERT_NO_THROW(GetValue(an_array_, Dimension{5}) = 9.9);
#endif
}

TEST_F(AnArrayForSTLContainerAccess, Returns0_5ForConstSTLContainerAccessAtFront) {  // NOLINT
  ASSERT_THAT(GetFront(an_array_), DoubleEq(0.5));
}

TEST_F(AnArrayForSTLContainerAccess, ThrowsExceptionForEmptyConstSTLContainerAccessAtFrontOnlyInDebugMode) {  // NOLINT
  std::array<int, 0> an_empty_array{};
#ifdef DEBUG
  ASSERT_THROW(GetFront(an_empty_array), std::out_of_range);
#else
  ASSERT_NO_THROW(GetFront(an_empty_array));
#endif
}

TEST_F(AnArrayForSTLContainerAccess, CanBeWritten10_1ToFrontViaSTLContainerAccess) {  // NOLINT
  GetFront(an_array_) = 10.1;
  ASSERT_THAT(GetValue(an_array_, 0), DoubleEq(10.1));
}

TEST_F(AnArrayForSTLContainerAccess,  // NOLINT
       ThrowsExceptionForWritingToEmptyFrontViaSTLContainerAccessOnlyInDebugMode) {
  std::array<int, 0> an_empty_array{};
#ifdef DEBUG
  ASSERT_THROW(GetFront(an_empty_array) = 9, std::out_of_range);
#else
  ASSERT_EXIT(GetFront(an_empty_array) = 9, testing::KilledBySignal(SIGSEGV), ".*");
#endif
}

TEST_F(AnArrayForSTLContainerAccess, ReturnsMinus5ForConstSTLContainerAccessAtBack) {  // NOLINT
  ASSERT_THAT(GetBack(an_array_), DoubleEq(-5));
}

TEST_F(AnArrayForSTLContainerAccess, ThrowsExceptionForEmptyConstSTLContainerAccessAtBackOnlyInDebugMode) {  // NOLINT
  std::array<int, 0> an_empty_array{};
#ifdef DEBUG
  ASSERT_THROW(GetBack(an_empty_array), std::out_of_range);
#else
  ASSERT_NO_THROW(GetBack(an_empty_array));
#endif
}

TEST_F(AnArrayForSTLContainerAccess, CanBeWritten10_1ToBackViaSTLContainerAccess) {  // NOLINT
  GetBack(an_array_) = 10.1;
  ASSERT_THAT(GetValue(an_array_, 4), DoubleEq(10.1));
}

TEST_F(AnArrayForSTLContainerAccess,  // NOLINT
       ThrowsExceptionForWritingToEmptyBackViaSTLContainerAccessOnlyInDebugMode) {
  std::array<int, 0> an_empty_array{};
#ifdef DEBUG
  ASSERT_THROW(GetBack(an_empty_array) = 9, std::out_of_range);
#else
  ASSERT_EXIT(GetBack(an_empty_array) = 9, testing::KilledBySignal(SIGSEGV), ".*");
#endif
}
