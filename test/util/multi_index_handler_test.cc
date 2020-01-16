/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/util/multi_index_handler.h"

using testing::Eq;
using testing::Test;

using namespace splinelib::src;

class AMultiIndexHandler1D : public Test {
 public:
  AMultiIndexHandler1D() {
    std::array<int, 1> number_of_knots = {10};
    multi_index_handler_1D_ = std::make_unique<util::MultiIndexHandler<1>>(number_of_knots);
  }

 protected:
  std::unique_ptr<util::MultiIndexHandler<1>> multi_index_handler_1D_;
};

TEST_F(AMultiIndexHandler1D, Returns1DIndex0AfterConstruction) {  // NOLINT
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(0));
}

TEST_F(AMultiIndexHandler1D, CanBeCopied) {  // NOLINT
  util::MultiIndexHandler<1> copied_multi_index_handler_1D = *multi_index_handler_1D_;
  ASSERT_THAT(copied_multi_index_handler_1D, Eq(*multi_index_handler_1D_));
}

TEST_F(AMultiIndexHandler1D, CanBeCreatedWithMoveConstructor) {  // NOLINT
  util::MultiIndexHandler<1> moved_multi_index_handler_1D(std::move(*multi_index_handler_1D_));
  ASSERT_THAT(moved_multi_index_handler_1D.GetNumberOfTotalMultiIndices(), Eq(10));
}

TEST_F(AMultiIndexHandler1D, CanBeAssigned) {  // NOLINT
  util::MultiIndexHandler<1> multi_index_handler_1D_to_assign{};
  multi_index_handler_1D_to_assign = *multi_index_handler_1D_;
  ASSERT_THAT(multi_index_handler_1D_to_assign, Eq(*multi_index_handler_1D_));
}

TEST_F(AMultiIndexHandler1D, CanBeMoveAssigned) {  // NOLINT
  util::MultiIndexHandler<1> multi_index_handler_1D_to_move_assign{};
  multi_index_handler_1D_to_move_assign = std::move(*multi_index_handler_1D_);
  ASSERT_THAT(multi_index_handler_1D_to_move_assign.GetNumberOfTotalMultiIndices(), Eq(10));
}

TEST_F(AMultiIndexHandler1D, Returns1DIndex1AfterUsingPostIncrementOperator) {  // NOLINT
  ASSERT_THAT(((*multi_index_handler_1D_)++).GetCurrent1DIndex(), Eq(0));
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(1));
}

TEST_F(AMultiIndexHandler1D, Returns1DIndex1AfterUsingPreIncrementOperator) {  // NOLINT
  ASSERT_THAT((++(*multi_index_handler_1D_)).GetCurrent1DIndex(), Eq(1));
}

TEST_F(AMultiIndexHandler1D, Returns1DIndex9AfterUsingPostDecrementOperator) {  // NOLINT
  ASSERT_THAT(((*multi_index_handler_1D_)--).GetCurrent1DIndex(), Eq(0));
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(9));
}

TEST_F(AMultiIndexHandler1D, Returns1DIndex9AfterUsingPreDecrementOperator) {  // NOLINT
  ASSERT_THAT((--(*multi_index_handler_1D_)).GetCurrent1DIndex(), Eq(9));
}

TEST_F(AMultiIndexHandler1D, Returns1DIndex1AfterAdding11) {  // NOLINT
  *multi_index_handler_1D_ + 11;
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(1));
}

TEST_F(AMultiIndexHandler1D, Returns1DIndex8AfterSubtracting2) {  // NOLINT
  *multi_index_handler_1D_ - 2;
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(8));
}

TEST_F(AMultiIndexHandler1D, CanBeIterated) {  // NOLINT
  int counter = 0;
  for (auto multi_index_handler_with_current_indices = multi_index_handler_1D_->begin();
       multi_index_handler_with_current_indices != multi_index_handler_1D_->end();
       ++multi_index_handler_with_current_indices, ++counter) {
    ASSERT_THAT(multi_index_handler_with_current_indices.GetCurrent1DIndex(), Eq(counter));
  }
  ASSERT_THAT(counter, Eq(multi_index_handler_1D_->GetNumberOfTotalMultiIndices()));
}

TEST_F(AMultiIndexHandler1D, ReturnsIndex5AtDimension0AfterSettingCurrentIndexTo1DIndex5) {  // NOLINT
  std::array<int, 1> current_multi_index = {5};
  multi_index_handler_1D_->SetCurrentIndex(current_multi_index);
  ASSERT_THAT((*multi_index_handler_1D_)[Dimension{0}], Eq(5));
}

TEST_F(AMultiIndexHandler1D, ThrowsInvalidArgumentWhenSettingCurrentIndexToMinus1) {  // NOLINT
  std::array<int, 1> current_multi_index = {-1};
  ASSERT_THROW(multi_index_handler_1D_->SetCurrentIndex(current_multi_index), std::invalid_argument);
}

TEST_F(AMultiIndexHandler1D, ThrowsInvalidArgumentWhenSettingCurrentIndexTo10) {  // NOLINT
  std::array<int, 1> current_multi_index = {10};
  ASSERT_THROW(multi_index_handler_1D_->SetCurrentIndex(current_multi_index), std::invalid_argument);
}

TEST_F(AMultiIndexHandler1D, Returns1DIndex5AfterSettingCurrentIndexTo5) {  // NOLINT
  multi_index_handler_1D_->SetCurrentIndex(5);
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(5));
}

TEST_F(AMultiIndexHandler1D, ThrowsInvalidArgumentWhenSettingCurrentIndexTo1DIndexMinus1) {  // NOLINT
  ASSERT_THROW(multi_index_handler_1D_->SetCurrentIndex(-1), std::invalid_argument);
}

TEST_F(AMultiIndexHandler1D, ThrowsInvalidArgumentWhenSettingCurrentIndexTo1DIndex10) {  // NOLINT
  ASSERT_THROW(multi_index_handler_1D_->SetCurrentIndex(10), std::invalid_argument);
}

TEST_F(AMultiIndexHandler1D, Returns1DIndex3AfterSettingCurrentIndexForDimension0To3) {  // NOLINT
  multi_index_handler_1D_->SetCurrentIndexForDimension(3, Dimension{0});
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(3));
}

TEST_F(AMultiIndexHandler1D, ThrowsInvalidArgumentWhenSettingCurrentIndexForTooSmallDimensionMinus1) {  // NOLINT
  ASSERT_THROW(multi_index_handler_1D_->SetCurrentIndexForDimension(3, Dimension{-1}), std::invalid_argument);
}

TEST_F(AMultiIndexHandler1D, ThrowsInvalidArgumentWhenSettingCurrentIndexForTooBigDimension1) {  // NOLINT
  ASSERT_THROW(multi_index_handler_1D_->SetCurrentIndexForDimension(3, Dimension{1}), std::invalid_argument);
}

TEST_F(AMultiIndexHandler1D, ThrowsInvalidArgumentWhenSettingCurrentIndexToTooSmallValueOfMinus1ForDimension0) {  // NOLINT
  ASSERT_THROW(multi_index_handler_1D_->SetCurrentIndexForDimension(-1, Dimension{0}), std::invalid_argument);
}

TEST_F(AMultiIndexHandler1D, ThrowsInvalidArgumentWhenSettingCurrentIndexToTooBigValueOf10ForDimension0) {  // NOLINT
  ASSERT_THROW(multi_index_handler_1D_->SetCurrentIndexForDimension(10, Dimension{0}), std::invalid_argument);
}

TEST_F(AMultiIndexHandler1D, Returns1DIndex1ForMultiIndex1WithLength2) {  // NOLINT
  std::array<int, 1> length_to_evaluate = {2};
  std::array<int, 1> multi_index_to_evaluate = {1};
  ASSERT_THAT(util::MultiIndexHandler<1>::Get1DIndex(length_to_evaluate, multi_index_to_evaluate), Eq(1));
}

TEST_F(AMultiIndexHandler1D, Returns1DIndex1ForMultiIndex1ButCurrent1DIndexIsStill0) {  // NOLINT
  std::array<int, 1> multi_index_to_evaluate = {1};
  ASSERT_THAT(multi_index_handler_1D_->Get1DIndex(multi_index_to_evaluate), Eq(1));
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(0));
}

TEST_F(AMultiIndexHandler1D, ReturnsComplementaryIndex4AfterSettingCurrentIndexTo5) {  // NOLINT
  multi_index_handler_1D_->SetCurrentIndex(5);
  ASSERT_THAT(multi_index_handler_1D_->GetComplementaryIndex()[0], Eq(4));
}

TEST_F(AMultiIndexHandler1D, Returns10AsNumberOfTotalMultiIndices) {  // NOLINT
  ASSERT_THAT(multi_index_handler_1D_->GetNumberOfTotalMultiIndices(), Eq(10));
}

TEST_F(AMultiIndexHandler1D, Returns0AfterCollapsingDimension0WithCurrentMultiIndex4) {  // NOLINT
  std::array<int, 1> current_multi_index = {4};
  multi_index_handler_1D_->SetCurrentIndex(current_multi_index);
  ASSERT_THAT(multi_index_handler_1D_->CollapseDimension(Dimension{0}), Eq(0));
}

TEST_F(AMultiIndexHandler1D, Returns1ForNumberOfMultiIndicesForCollapsedDimension0) {  // NOLINT
  ASSERT_THAT(multi_index_handler_1D_->GetNumberOfMultiIndicesForCollapsedDimension(Dimension{0}), Eq(1));
}

TEST_F(AMultiIndexHandler1D, Returns1ForCurrentSliceSizeUpToDimension0) {  // NOLINT
  ASSERT_THAT(multi_index_handler_1D_->GetCurrentSliceSize(Dimension{0}), Eq(1));
}

TEST_F(AMultiIndexHandler1D, Returns1ForCurrentSliceComplementFillUpToDimension0WithCurrentMultiIndex0) {  // NOLINT
  ASSERT_THAT(multi_index_handler_1D_->GetCurrentSliceComplementFill(Dimension{0}), Eq(1));
}

class AMultiIndexHandler2D : public Test {
 public:
  AMultiIndexHandler2D() {
    std::array<int, 2> number_of_knots_per_dimension = {9, 3};
    multi_index_handler_2D_ = std::make_unique<util::MultiIndexHandler<2>>(number_of_knots_per_dimension);
  }

 protected:
  std::unique_ptr<util::MultiIndexHandler<2>> multi_index_handler_2D_;
};

TEST_F(AMultiIndexHandler2D, Returns1DIndex0AfterConstruction) {  // NOLINT
  ASSERT_THAT(multi_index_handler_2D_->GetCurrent1DIndex(), Eq(0));
}

TEST_F(AMultiIndexHandler2D, CanBeCopied) {  // NOLINT
  util::MultiIndexHandler<2> copied_multi_index_handler_2D = *multi_index_handler_2D_;
  ASSERT_THAT(copied_multi_index_handler_2D, Eq(*multi_index_handler_2D_));
}

TEST_F(AMultiIndexHandler2D, CanBeAssigned) {  // NOLINT
  util::MultiIndexHandler<2> multi_index_handler_2D_to_be_assigned{};
  multi_index_handler_2D_to_be_assigned = *multi_index_handler_2D_;
  ASSERT_THAT(multi_index_handler_2D_to_be_assigned, Eq(*multi_index_handler_2D_));
}

TEST_F(AMultiIndexHandler2D, ReturnsIndex1And0AfterUsingPostIncrementOperator) {  // NOLINT
  ASSERT_THAT(((*multi_index_handler_2D_)++).GetCurrent1DIndex(), Eq(0));
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{0}], Eq(1));
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{1}], Eq(0));
}

TEST_F(AMultiIndexHandler2D, Returns1DIndex1AfterUsingPreIncrementOperator) {  // NOLINT
  ASSERT_THAT((++(*multi_index_handler_2D_)).GetCurrent1DIndex(), Eq(1));
}

TEST_F(AMultiIndexHandler2D, ReturnsIndex8And2AfterUsingPostDecrementOperator) {  // NOLINT
  ASSERT_THAT(((*multi_index_handler_2D_)--).GetCurrent1DIndex(), Eq(0));
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{0}], Eq(8));
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{1}], Eq(2));
}

TEST_F(AMultiIndexHandler2D, Returns1DIndex26AfterUsingPreDecrementOperator) {  // NOLINT
  ASSERT_THAT((--(*multi_index_handler_2D_)).GetCurrent1DIndex(), Eq(26));
}

TEST_F(AMultiIndexHandler2D, ReturnsIndex1And1AfterAdding10) {  // NOLINT
  *multi_index_handler_2D_ + 10;
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{0}], Eq(1));
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{1}], Eq(1));
}

TEST_F(AMultiIndexHandler2D, ReturnsIndex6And1AfterSubtracting12) {  // NOLINT
  *multi_index_handler_2D_ - 12;
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{0}], Eq(6));
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{1}], Eq(1));
}

TEST_F(AMultiIndexHandler2D, ReturnsIndex5And0AfterSubtractionFromIndex5And1WithConstantDimension0) {  // NOLINT
  multi_index_handler_2D_->SetCurrentIndex({5, 1});
  multi_index_handler_2D_->SubtractWithConstantDimension(Dimension{0});
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{0}], Eq(5));
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{1}], Eq(0));
}

TEST_F(AMultiIndexHandler2D, ReturnsIndex4And1AfterSubtractionFromIndex5And1WithConstantDimension1) {  // NOLINT
  multi_index_handler_2D_->SetCurrentIndex({5, 1});
  multi_index_handler_2D_->SubtractWithConstantDimension(Dimension{1});
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{0}], Eq(4));
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{1}], Eq(1));
}

TEST_F(AMultiIndexHandler2D, CanBeIterated) {  // NOLINT
  int counter = 0;
  for (auto multi_index_handler_with_current_indices = multi_index_handler_2D_->begin();
       multi_index_handler_with_current_indices != multi_index_handler_2D_->end();
       ++multi_index_handler_with_current_indices, ++counter) {
    ASSERT_THAT(multi_index_handler_with_current_indices.GetCurrent1DIndex(), Eq(counter));
  }
  ASSERT_THAT(counter, Eq(multi_index_handler_2D_->GetNumberOfTotalMultiIndices()));
}

TEST_F(AMultiIndexHandler2D, Returns1DIndex21AfterSettingCurrentIndexToIndex3And2) {  // NOLINT
  std::array<int, 2> current_multi_index = {3, 2};
  multi_index_handler_2D_->SetCurrentIndex(current_multi_index);
  ASSERT_THAT(multi_index_handler_2D_->GetCurrent1DIndex(), Eq(21));
}

TEST_F(AMultiIndexHandler2D, ThrowsInvalidArgumentWhenSettingCurrentIndexTo3AndMinus1) {  // NOLINT
  std::array<int, 2> current_multi_index = {3, -1};
  ASSERT_THROW(multi_index_handler_2D_->SetCurrentIndex(current_multi_index), std::invalid_argument);
}

TEST_F(AMultiIndexHandler2D, ThrowsInvalidArgumentWhenSettingCurrentIndexTo9And2) {  // NOLINT
  std::array<int, 2> current_multi_index = {9, 2};
  ASSERT_THROW(multi_index_handler_2D_->SetCurrentIndex(current_multi_index), std::invalid_argument);
}

TEST_F(AMultiIndexHandler2D, ReturnsIndex3And2AfterSettingCurrentIndexTo1DIndex21) {  // NOLINT
  multi_index_handler_2D_->SetCurrentIndex(21);
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{0}], Eq(3));
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{1}], Eq(2));
}

TEST_F(AMultiIndexHandler2D, ThrowsInvalidArgumentWhenSettingCurrentIndexTo1DIndexMinus1) {  // NOLINT
  ASSERT_THROW(multi_index_handler_2D_->SetCurrentIndex(-1), std::invalid_argument);
}

TEST_F(AMultiIndexHandler2D, ThrowsInvalidArgumentWhenSettingCurrentIndexTo1DIndex27) {  // NOLINT
  ASSERT_THROW(multi_index_handler_2D_->SetCurrentIndex(27), std::invalid_argument);
}

TEST_F(AMultiIndexHandler2D, Returns1DIndex3AfterSettingCurrentIndexForDimension0To3) {  // NOLINT
  multi_index_handler_2D_->SetCurrentIndexForDimension(3, Dimension{0});
  ASSERT_THAT(multi_index_handler_2D_->GetCurrent1DIndex(), Eq(3));
}

TEST_F(AMultiIndexHandler2D, Returns1DIndex9AfterSettingCurrentIndexForDimension1To1) {  // NOLINT
  multi_index_handler_2D_->SetCurrentIndexForDimension(1, Dimension{1});
  ASSERT_THAT(multi_index_handler_2D_->GetCurrent1DIndex(), Eq(9));
}

TEST_F(AMultiIndexHandler2D, ThrowsInvalidArgumentWhenSettingCurrentIndexForTooSmallDimensionMinus1) {  // NOLINT
  ASSERT_THROW(multi_index_handler_2D_->SetCurrentIndexForDimension(3, Dimension{-1}), std::invalid_argument);
}

TEST_F(AMultiIndexHandler2D, ThrowsInvalidArgumentWhenSettingCurrentIndexForTooBigDimension2) {  // NOLINT
  ASSERT_THROW(multi_index_handler_2D_->SetCurrentIndexForDimension(3, Dimension{2}), std::invalid_argument);
}

TEST_F(AMultiIndexHandler2D, ThrowsInvalidArgumentWhenSettingCurrentIndexToTooSmallValueOfMinus1ForDimension0And1) {  // NOLINT
  ASSERT_THROW(multi_index_handler_2D_->SetCurrentIndexForDimension(-1, Dimension{0}), std::invalid_argument);
  ASSERT_THROW(multi_index_handler_2D_->SetCurrentIndexForDimension(-1, Dimension{1}), std::invalid_argument);
}

TEST_F(AMultiIndexHandler2D, ThrowsInvalidArgumentWhenSettingCurrentIndexToTooBigValueOf9And3ForDimension0And1) {  // NOLINT
  ASSERT_THROW(multi_index_handler_2D_->SetCurrentIndexForDimension(9, Dimension{0}), std::invalid_argument);
  ASSERT_THROW(multi_index_handler_2D_->SetCurrentIndexForDimension(3, Dimension{1}), std::invalid_argument);
}

TEST_F(AMultiIndexHandler2D, Returns1DIndex6ForMultiIndex2And1WithLength4And3) {  // NOLINT
  std::array<int, 2> length_to_evaluate = {4, 3};
  std::array<int, 2> multi_index_to_evaluate = {2, 1};
  ASSERT_THAT(util::MultiIndexHandler<2>::Get1DIndex(length_to_evaluate, multi_index_to_evaluate), Eq(6));
}

TEST_F(AMultiIndexHandler2D, Returns1DIndex11ForMultiIndex2And1ButCurrent1DIndexIsStill0) {  // NOLINT
  std::array<int, 2> multi_index_to_evaluate = {2, 1};
  ASSERT_THAT(multi_index_handler_2D_->Get1DIndex(multi_index_to_evaluate), Eq(11));
  ASSERT_THAT(multi_index_handler_2D_->GetCurrent1DIndex(), Eq(0));
}

TEST_F(AMultiIndexHandler2D, ReturnsComplementaryIndex3And1AfterSettingCurrentIndexTo5And1) {  // NOLINT
  std::array<int, 2> current_multi_index = {5, 1};
  multi_index_handler_2D_->SetCurrentIndex(current_multi_index);
  ASSERT_THAT(multi_index_handler_2D_->GetComplementaryIndex()[0], Eq(3));
  ASSERT_THAT(multi_index_handler_2D_->GetComplementaryIndex()[1], Eq(1));
}

TEST_F(AMultiIndexHandler2D, Returns30AsNumberOfTotalMultiIndices) {  // NOLINT
  ASSERT_THAT(multi_index_handler_2D_->GetNumberOfTotalMultiIndices(), Eq(27));
}

TEST_F(AMultiIndexHandler2D, Returns1DIndex1AfterCollapsingDimension0WithCurrentMultiIndex4And1) {  // NOLINT
  std::array<int, 2> current_multi_index = {4, 1};
  multi_index_handler_2D_->SetCurrentIndex(current_multi_index);
  ASSERT_THAT(multi_index_handler_2D_->CollapseDimension(Dimension{0}), Eq(1));
}

TEST_F(AMultiIndexHandler2D, Returns1DIndex4AfterCollapsingDimension1WithCurrentMultiIndex4And1) {  // NOLINT
  std::array<int, 2> current_multi_index = {4, 1};
  multi_index_handler_2D_->SetCurrentIndex(current_multi_index);
  ASSERT_THAT(multi_index_handler_2D_->CollapseDimension(Dimension{1}), Eq(4));
}

TEST_F(AMultiIndexHandler2D, Returns3ForNumberOfMultiIndicesForCollapsedDimension0) {  // NOLINT
  ASSERT_THAT(multi_index_handler_2D_->GetNumberOfMultiIndicesForCollapsedDimension(Dimension{0}), Eq(3));
}

TEST_F(AMultiIndexHandler2D, Returns9ForNumberOfMultiIndicesForCollapsedDimension1) {  // NOLINT
  ASSERT_THAT(multi_index_handler_2D_->GetNumberOfMultiIndicesForCollapsedDimension(Dimension{1}), Eq(9));
}

TEST_F(AMultiIndexHandler2D, Returns1ForCurrentSliceSizeUpToDimension0) {  // NOLINT
  ASSERT_THAT(multi_index_handler_2D_->GetCurrentSliceSize(Dimension{0}), Eq(1));
}

TEST_F(AMultiIndexHandler2D, Returns9ForCurrentSliceSizeUpToDimension1) {  // NOLINT
  ASSERT_THAT(multi_index_handler_2D_->GetCurrentSliceSize(Dimension{1}), Eq(9));
}

TEST_F(AMultiIndexHandler2D, Returns1ForCurrentSliceComplementFillUpToDimension0WithCurrentMultiIndex0And0) {  // NOLINT
  ASSERT_THAT(multi_index_handler_2D_->GetCurrentSliceComplementFill(Dimension{0}), Eq(1));
}

TEST_F(AMultiIndexHandler2D, Returns8ForCurrentSliceComplementFillUpToDimension1WithCurrentMultiIndex1And0) {  // NOLINT
  multi_index_handler_2D_->SetCurrentIndex({1, 0});
  ASSERT_THAT(multi_index_handler_2D_->GetCurrentSliceComplementFill(Dimension{1}), Eq(8));
}

class AMultiIndexHandler3D : public Test {
 public:
  AMultiIndexHandler3D() {
    std::array<int, 3> number_of_knots_per_dimension = {4, 3, 5};
    multi_index_handler_3D_ = std::make_unique<util::MultiIndexHandler<3>>(number_of_knots_per_dimension);
  }
 protected:
  std::unique_ptr<util::MultiIndexHandler<3>> multi_index_handler_3D_;
};

TEST_F(AMultiIndexHandler3D, Returns1DIndex0AfterConstruction) {  // NOLINT
  ASSERT_THAT(multi_index_handler_3D_->GetCurrent1DIndex(), Eq(0));
}

TEST_F(AMultiIndexHandler3D, CanBeCopied) {  // NOLINT
  util::MultiIndexHandler<3> copied_multi_index_handler_3D = *multi_index_handler_3D_;
  ASSERT_THAT(copied_multi_index_handler_3D, Eq(*multi_index_handler_3D_));
}

TEST_F(AMultiIndexHandler3D, CanBeAssigned) {  // NOLINT
  util::MultiIndexHandler<3> multi_index_handler_3D_to_be_assigned{};
  multi_index_handler_3D_to_be_assigned = *multi_index_handler_3D_;
  ASSERT_THAT(multi_index_handler_3D_to_be_assigned, Eq(*multi_index_handler_3D_));
}

TEST_F(AMultiIndexHandler3D, ReturnsIndex1And0And0AfterUsingPostIncrementOperator) {  // NOLINT
  ASSERT_THAT(((*multi_index_handler_3D_)++).GetCurrent1DIndex(), Eq(0));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{0}], Eq(1));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{1}], Eq(0));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{2}], Eq(0));
}

TEST_F(AMultiIndexHandler3D, Returns1DIndex1AfterUsingPreIncrementOperator) {  // NOLINT
  ASSERT_THAT((++(*multi_index_handler_3D_)).GetCurrent1DIndex(), Eq(1));
}

TEST_F(AMultiIndexHandler3D, ReturnsIndex3And2And4AfterUsingPostDecrementOperator) {  // NOLINT
  ASSERT_THAT(((*multi_index_handler_3D_)--).GetCurrent1DIndex(), Eq(0));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{0}], Eq(3));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{1}], Eq(2));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{2}], Eq(4));
}

TEST_F(AMultiIndexHandler3D, Returns1DIndex59AfterUsingPreDecrementOperator) {  // NOLINT
  ASSERT_THAT((--(*multi_index_handler_3D_)).GetCurrent1DIndex(), Eq(59));
}

TEST_F(AMultiIndexHandler3D, ReturnsIndex3And0And1AfterAdding15) {  // NOLINT
  *multi_index_handler_3D_ + 15;
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{0}], Eq(3));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{1}], Eq(0));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{2}], Eq(1));
}

TEST_F(AMultiIndexHandler3D, ReturnsIndex6And1AfterSubtracting22) {  // NOLINT
  *multi_index_handler_3D_ - 22;
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{0}], Eq(2));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{1}], Eq(0));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{2}], Eq(3));
}

TEST_F(AMultiIndexHandler3D, ReturnsIndex2And2And3AfterSubtractionFromIndex2And0And4WithConstantDimension0) {  // NOLINT
  multi_index_handler_3D_->SetCurrentIndex({2, 0, 4});
  multi_index_handler_3D_->SubtractWithConstantDimension(Dimension{0});
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{0}], Eq(2));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{1}], Eq(2));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{2}], Eq(3));
}

TEST_F(AMultiIndexHandler3D, ReturnsIndex1And0And4AfterSubtractionFromIndex2And0And4WithConstantDimension1) {  // NOLINT
  multi_index_handler_3D_->SetCurrentIndex({2, 0, 4});
  multi_index_handler_3D_->SubtractWithConstantDimension(Dimension{1});
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{0}], Eq(1));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{1}], Eq(0));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{2}], Eq(4));
}

TEST_F(AMultiIndexHandler3D, ReturnsIndex1And0And4AfterSubtractionFromIndex2And0And4WithConstantDimension2) {  // NOLINT
  multi_index_handler_3D_->SetCurrentIndex({2, 0, 4});
  multi_index_handler_3D_->SubtractWithConstantDimension(Dimension{2});
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{0}], Eq(1));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{1}], Eq(0));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{2}], Eq(4));
}

TEST_F(AMultiIndexHandler3D, CanBeIterated) {  // NOLINT
  int counter = 0;
  for (auto multi_index_handler_with_current_indices = multi_index_handler_3D_->begin();
       multi_index_handler_with_current_indices != multi_index_handler_3D_->end();
       ++multi_index_handler_with_current_indices, ++counter) {
    ASSERT_THAT(multi_index_handler_with_current_indices.GetCurrent1DIndex(), Eq(counter));
  }
  ASSERT_THAT(counter, Eq(multi_index_handler_3D_->GetNumberOfTotalMultiIndices()));
}

TEST_F(AMultiIndexHandler3D, Returns1DIndex29AfterSettingCurrentIndexToIndex1And1And2) {  // NOLINT
  std::array<int, 3> current_multi_index = {1, 1, 2};
  multi_index_handler_3D_->SetCurrentIndex(current_multi_index);
  ASSERT_THAT(multi_index_handler_3D_->GetCurrent1DIndex(), Eq(29));
}

TEST_F(AMultiIndexHandler3D, ThrowsInvalidArgumentWhenSettingCurrentIndexTo3AndMinus1And3) {  // NOLINT
  std::array<int, 3> current_multi_index = {3, -1, 3};
  ASSERT_THROW(multi_index_handler_3D_->SetCurrentIndex(current_multi_index), std::invalid_argument);
}

TEST_F(AMultiIndexHandler3D, ThrowsInvalidArgumentWhenSettingCurrentIndexTo3And3And2) {  // NOLINT
  std::array<int, 3> current_multi_index = {3, 3, 2};
  ASSERT_THROW(multi_index_handler_3D_->SetCurrentIndex(current_multi_index), std::invalid_argument);
}

TEST_F(AMultiIndexHandler3D, ReturnsIndex1And1And2AfterSettingCurrentIndexTo1DIndex29) {  // NOLINT
  multi_index_handler_3D_->SetCurrentIndex(29);
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{0}], Eq(1));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{1}], Eq(1));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{2}], Eq(2));
}

TEST_F(AMultiIndexHandler3D, ThrowsInvalidArgumentWhenSettingCurrentIndexTo1DIndexMinus1) {  // NOLINT
  ASSERT_THROW(multi_index_handler_3D_->SetCurrentIndex(-1), std::invalid_argument);
}

TEST_F(AMultiIndexHandler3D, ThrowsInvalidArgumentWhenSettingCurrentIndexTo1DIndex60) {  // NOLINT
  ASSERT_THROW(multi_index_handler_3D_->SetCurrentIndex(60), std::invalid_argument);
}

TEST_F(AMultiIndexHandler3D, Returns1DIndex3AfterSettingCurrentIndexForDimension0To3) {  // NOLINT
  multi_index_handler_3D_->SetCurrentIndexForDimension(3, Dimension{0});
  ASSERT_THAT(multi_index_handler_3D_->GetCurrent1DIndex(), Eq(3));
}

TEST_F(AMultiIndexHandler3D, Returns1DIndex4AfterSettingCurrentIndexForDimension1To1) {  // NOLINT
  multi_index_handler_3D_->SetCurrentIndexForDimension(1, Dimension{1});
  ASSERT_THAT(multi_index_handler_3D_->GetCurrent1DIndex(), Eq(4));
}

TEST_F(AMultiIndexHandler3D, Returns1DIndex24AfterSettingCurrentIndexForDimension2To2) {  // NOLINT
  multi_index_handler_3D_->SetCurrentIndexForDimension(2, Dimension{2});
  ASSERT_THAT(multi_index_handler_3D_->GetCurrent1DIndex(), Eq(24));
}

TEST_F(AMultiIndexHandler3D, ThrowsInvalidArgumentWhenSettingCurrentIndexForTooSmallDimensionMinus1) {  // NOLINT
  ASSERT_THROW(multi_index_handler_3D_->SetCurrentIndexForDimension(3, Dimension{-1}), std::invalid_argument);
}

TEST_F(AMultiIndexHandler3D, ThrowsInvalidArgumentWhenSettingCurrentIndexForTooBigDimension3) {  // NOLINT
  ASSERT_THROW(multi_index_handler_3D_->SetCurrentIndexForDimension(3, Dimension{3}), std::invalid_argument);
}

TEST_F(AMultiIndexHandler3D, ThrowsInvalidArgumentWhenSettingCurrentIndexToTooSmallValueMinus1ForDimension0And1And2) {  // NOLINT
  ASSERT_THROW(multi_index_handler_3D_->SetCurrentIndexForDimension(-1, Dimension{0}), std::invalid_argument);
  ASSERT_THROW(multi_index_handler_3D_->SetCurrentIndexForDimension(-1, Dimension{1}), std::invalid_argument);
  ASSERT_THROW(multi_index_handler_3D_->SetCurrentIndexForDimension(-1, Dimension{2}), std::invalid_argument);
}

TEST_F(AMultiIndexHandler3D, ThrowsInvalidArgumentWhenSettingCurrentIndexToTooBigValue4And3And5ForDimension0And1And2) {  // NOLINT
  ASSERT_THROW(multi_index_handler_3D_->SetCurrentIndexForDimension(4, Dimension{0}), std::invalid_argument);
  ASSERT_THROW(multi_index_handler_3D_->SetCurrentIndexForDimension(3, Dimension{1}), std::invalid_argument);
  ASSERT_THROW(multi_index_handler_3D_->SetCurrentIndexForDimension(5, Dimension{2}), std::invalid_argument);
}

TEST_F(AMultiIndexHandler3D, Returns1DIndex17ForMultiIndex2And1And1WithLength3And4And2) {  // NOLINT
  std::array<int, 3> length_to_evaluate = {3, 4, 2};
  std::array<int, 3> multi_index_to_evaluate = {2, 1, 1};
  ASSERT_THAT(util::MultiIndexHandler<3>::Get1DIndex(length_to_evaluate, multi_index_to_evaluate), Eq(17));
}

TEST_F(AMultiIndexHandler3D, Returns1DIndex18ForMultiIndex2And1And1ButCurrent1DIndexIsStill0) {  // NOLINT
  std::array<int, 3> multi_index_to_evaluate = {2, 1, 1};
  ASSERT_THAT(multi_index_handler_3D_->Get1DIndex(multi_index_to_evaluate), Eq(18));
  ASSERT_THAT(multi_index_handler_3D_->GetCurrent1DIndex(), Eq(0));
}

TEST_F(AMultiIndexHandler3D, ReturnsComplementaryIndex1And1And3AfterSettingCurrentIndexTo2And1And1) {  // NOLINT
  std::array<int, 3> current_multi_index = {2, 1, 1};
  multi_index_handler_3D_->SetCurrentIndex(current_multi_index);
  ASSERT_THAT(multi_index_handler_3D_->GetComplementaryIndex()[0], Eq(1));
  ASSERT_THAT(multi_index_handler_3D_->GetComplementaryIndex()[1], Eq(1));
  ASSERT_THAT(multi_index_handler_3D_->GetComplementaryIndex()[2], Eq(3));
}

TEST_F(AMultiIndexHandler3D, Returns60AsNumberOfTotalMultiIndices) {  // NOLINT
  ASSERT_THAT(multi_index_handler_3D_->GetNumberOfTotalMultiIndices(), Eq(60));
}

TEST_F(AMultiIndexHandler3D, Returns1DIndex7AfterCollapsingDimension0WithCurrentMultiIndex3And1And2) {  // NOLINT
  std::array<int, 3> current_multi_index = {3, 1, 2};
  multi_index_handler_3D_->SetCurrentIndex(current_multi_index);
  ASSERT_THAT(multi_index_handler_3D_->CollapseDimension(Dimension{0}), Eq(7));
}

TEST_F(AMultiIndexHandler3D, Returns1DIndex11AfterCollapsingDimension1WithCurrentMultiIndex3And1And2) {  // NOLINT
  std::array<int, 3> current_multi_index = {3, 1, 2};
  multi_index_handler_3D_->SetCurrentIndex(current_multi_index);
  ASSERT_THAT(multi_index_handler_3D_->CollapseDimension(Dimension{1}), Eq(11));
}

TEST_F(AMultiIndexHandler3D, Returns1DIndex7AfterCollapsingDimension2WithCurrentMultiIndex3And1And2) {  // NOLINT
  std::array<int, 3> current_multi_index = {3, 1, 2};
  multi_index_handler_3D_->SetCurrentIndex(current_multi_index);
  ASSERT_THAT(multi_index_handler_3D_->CollapseDimension(Dimension{2}), Eq(7));
}

TEST_F(AMultiIndexHandler3D, Returns15ForNumberOfMultiIndicesForCollapsedDimension0) {  // NOLINT
  ASSERT_THAT(multi_index_handler_3D_->GetNumberOfMultiIndicesForCollapsedDimension(Dimension{0}), Eq(15));
}

TEST_F(AMultiIndexHandler3D, Returns20ForNumberOfMultiIndicesForCollapsedDimension1) {  // NOLINT
  ASSERT_THAT(multi_index_handler_3D_->GetNumberOfMultiIndicesForCollapsedDimension(Dimension{1}), Eq(20));
}

TEST_F(AMultiIndexHandler3D, Returns12ForNumberOfMultiIndicesForCollapsedDimension2) {  // NOLINT
  ASSERT_THAT(multi_index_handler_3D_->GetNumberOfMultiIndicesForCollapsedDimension(Dimension{2}), Eq(12));
}

TEST_F(AMultiIndexHandler3D, Returns1ForCurrentSliceSizeUpToDimension0) {  // NOLINT
  ASSERT_THAT(multi_index_handler_3D_->GetCurrentSliceSize(Dimension{0}), Eq(1));
}

TEST_F(AMultiIndexHandler3D, Returns4ForCurrentSliceSizeUpToDimension1) {  // NOLINT
  ASSERT_THAT(multi_index_handler_3D_->GetCurrentSliceSize(Dimension{1}), Eq(4));
}

TEST_F(AMultiIndexHandler3D, Returns12ForCurrentSliceSizeUpToDimension2) {  // NOLINT
  ASSERT_THAT(multi_index_handler_3D_->GetCurrentSliceSize(Dimension{2}), Eq(12));
}

TEST_F(AMultiIndexHandler3D, Returns1ForCurrentSliceComplementFillUpToDimension0) {  // NOLINT
  ASSERT_THAT(multi_index_handler_3D_->GetCurrentSliceComplementFill(Dimension{0}), Eq(1));
}

TEST_F(AMultiIndexHandler3D, Returns3ForCurrentSliceComplementFillUpToDimension1WithCurrentMultiIndex1And1And1) {  // NOLINT
  multi_index_handler_3D_->SetCurrentIndex({1, 1, 1});
  ASSERT_THAT(multi_index_handler_3D_->GetCurrentSliceComplementFill(Dimension{1}), Eq(3));
}

TEST_F(AMultiIndexHandler3D, Returns7ForCurrentSliceComplementFillUpToDimension2WithCurrentMultiIndex1And1And1) {  // NOLINT
  multi_index_handler_3D_->SetCurrentIndex({1, 1, 1});
  ASSERT_THAT(multi_index_handler_3D_->GetCurrentSliceComplementFill(Dimension{2}), Eq(7));
}
