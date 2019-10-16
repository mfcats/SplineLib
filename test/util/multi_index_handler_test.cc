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

using testing::Test;
using testing::Eq;

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

TEST_F(AMultiIndexHandler1D, CanBeCopied) { // NOLINT
  util::MultiIndexHandler<1> copied_multi_index_handler_1D = *multi_index_handler_1D_;
  ASSERT_THAT(copied_multi_index_handler_1D, Eq(*multi_index_handler_1D_));
}

TEST_F(AMultiIndexHandler1D, CanBeAssigned) { // NOLINT
  util::MultiIndexHandler<1> multi_index_handler_1D_to_be_assigned{};
  multi_index_handler_1D_to_be_assigned = *multi_index_handler_1D_;
  ASSERT_THAT(multi_index_handler_1D_to_be_assigned, Eq(*multi_index_handler_1D_));
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
  (*multi_index_handler_1D_) + 11;
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(1));
}

TEST_F(AMultiIndexHandler1D, Returns1DIndex8AfterSubtracting2) {  // NOLINT
  (*multi_index_handler_1D_) - 2;
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(8));
}

TEST_F(AMultiIndexHandler1D, ReturnsIndex5AtDimension0AfterSettingCurrentIndexTo1DIndex5) {  // NOLINT
  std::array<int, 1> current_1D_index = {5};
  multi_index_handler_1D_->SetCurrentIndex(current_1D_index);
  ASSERT_THAT((*multi_index_handler_1D_)[Dimension{0}], Eq(5));
}

TEST_F(AMultiIndexHandler1D, Returns1DIndex5AfterSettingCurrentIndexTo5) {  // NOLINT
  multi_index_handler_1D_->SetCurrentIndex(5);
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(5));
}

TEST_F(AMultiIndexHandler1D, Returns1DIndex1ForMultiIndex1WithLength2) {  // NOLINT
  std::array<int, 1> length_to_evaluate = {2};
  std::array<int, 1> multi_index_to_evaluate = {1};
  ASSERT_THAT(multi_index_handler_1D_->Get1DIndex(length_to_evaluate, multi_index_to_evaluate), Eq(1));
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

TEST_F(AMultiIndexHandler2D, CanBeCopied) { // NOLINT
  util::MultiIndexHandler<2> copied_multi_index_handler_2D = *multi_index_handler_2D_;
  ASSERT_THAT(copied_multi_index_handler_2D, Eq(*multi_index_handler_2D_));
}

TEST_F(AMultiIndexHandler2D, CanBeAssigned) { // NOLINT
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
  (*multi_index_handler_2D_) + 10;
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{0}], Eq(1));
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{1}], Eq(1));
}

TEST_F(AMultiIndexHandler2D, ReturnsIndex6And1AfterSubtracting12) {  // NOLINT
  (*multi_index_handler_2D_) - 12;
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{0}], Eq(6));
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{1}], Eq(1));
}

TEST_F(AMultiIndexHandler2D, Returns1DIndex21AfterSettingCurrentIndexToIndex3And2) {  // NOLINT
  std::array<int, 2> current_multi_index = {3, 2};
  multi_index_handler_2D_->SetCurrentIndex(current_multi_index);
  ASSERT_THAT(multi_index_handler_2D_->GetCurrent1DIndex(), Eq(21));
}

TEST_F(AMultiIndexHandler2D, ReturnsIndex3And2AfterSettingCurrentIndexTo1DIndex21) {  // NOLINT
  multi_index_handler_2D_->SetCurrentIndex(21);
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{0}], Eq(3));
  ASSERT_THAT((*multi_index_handler_2D_)[Dimension{1}], Eq(2));
}

TEST_F(AMultiIndexHandler2D, Returns1DIndex6ForMultiIndex2And1WithLength4And3) {  // NOLINT
  std::array<int, 2> length_to_evaluate = {4, 3};
  std::array<int, 2> multi_index_to_evaluate = {2, 1};
  ASSERT_THAT(multi_index_handler_2D_->Get1DIndex(length_to_evaluate, multi_index_to_evaluate), Eq(6));
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

class MultiIndexHandler3D : public Test {
 public:
  MultiIndexHandler3D() {
    std::array<int, 3> last_knot_offset_3D = {4, 3, 5};
    multi_index_handler_3D_ = std::make_unique<util::MultiIndexHandler<3>>(last_knot_offset_3D);
  }
 protected:
  std::unique_ptr<util::MultiIndexHandler<3>> multi_index_handler_3D_;
};

TEST_F(MultiIndexHandler3D, Returns1DIndex0AfterConstruction) { // NOLINT
  ASSERT_THAT(multi_index_handler_3D_->GetCurrent1DIndex(), Eq(0));
}

TEST_F(MultiIndexHandler3D, Returns1DIndex30AfterSettingCurrentIndexTo2And1And2) { // NOLINT
  std::array<int, 3> current_3d_index = {2, 1, 2};
  multi_index_handler_3D_->SetCurrentIndex(current_3d_index);
  ASSERT_THAT(multi_index_handler_3D_->GetCurrent1DIndex(), Eq(30));
}

TEST_F(MultiIndexHandler3D, ReturnsIndex3And1And0AfterAdding1DIndex7) { // NOLINT
  *multi_index_handler_3D_ + 7;
  ASSERT_THAT(multi_index_handler_3D_->GetCurrent1DIndex(), Eq(7));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{0}], Eq(3));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{1}], Eq(1));
  ASSERT_THAT((*multi_index_handler_3D_)[Dimension{2}], Eq(0));
}

TEST_F(MultiIndexHandler3D, ReturnsComplementaryIndex0And1And4AfterAdding1DIndex7) { // NOLINT
  *multi_index_handler_3D_ + 7;
  auto complementary_multi_index = multi_index_handler_3D_->GetComplementaryIndex();
  ASSERT_THAT(multi_index_handler_3D_->GetCurrent1DIndex(), Eq(7));
  ASSERT_THAT(complementary_multi_index[0], Eq(0));
  ASSERT_THAT(complementary_multi_index[1], Eq(1));
  ASSERT_THAT(complementary_multi_index[2], Eq(4));
}

TEST_F(MultiIndexHandler3D, ReturnsCorrectNumberOf60DifferentMultiIndices) { // NOLINT
  ASSERT_THAT(multi_index_handler_3D_->GetNumberOfTotalMultiIndices(), Eq(60));
}

TEST_F(MultiIndexHandler3D, ReturnsCorrect1DIndexOf14CollapsingDimension1) { // NOLINT
  std::array<int, 3> indices = {2, 1, 3};
  multi_index_handler_3D_->SetCurrentIndex(indices);
  ASSERT_THAT(multi_index_handler_3D_->CollapseDimension(Dimension{1}), Eq(14));
}
