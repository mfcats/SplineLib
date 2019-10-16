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
  util::MultiIndexHandler<1> assigned_multi_index_handler_1D{};
  assigned_multi_index_handler_1D = *multi_index_handler_1D_;
  ASSERT_THAT(assigned_multi_index_handler_1D, Eq(*multi_index_handler_1D_));
}

TEST_F(AMultiIndexHandler1D, ReturnsIndex1AfterUsingPostIncrementOperator) {  // NOLINT
  (*multi_index_handler_1D_)++;
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(1));
}

TEST_F(AMultiIndexHandler1D, ReturnsIndex1AfterUsingPreIncrementOperator) {  // NOLINT
  ++(*multi_index_handler_1D_);
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(1));
}

TEST_F(AMultiIndexHandler1D, ReturnsIndex9AfterUsingPostDecrementOperator) {  // NOLINT
  (*multi_index_handler_1D_)--;
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(9));
}

TEST_F(AMultiIndexHandler1D, ReturnsIndex1AfterUsingPreDecrementOperator) {  // NOLINT
  --(*multi_index_handler_1D_);
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(9));
}

TEST_F(AMultiIndexHandler1D, ReturnsIndex10AfterAdding10) {  // NOLINT
  (*multi_index_handler_1D_) = (*multi_index_handler_1D_) + 10;
  ASSERT_THAT(multi_index_handler_1D_->GetCurrent1DIndex(), Eq(0));
}

TEST_F(AMultiIndexHandler1D, ReturnsIndex8AfterSubtracting2) {  // NOLINT
  (*multi_index_handler_1D_) = (*multi_index_handler_1D_) - 2;
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
  std::array<int, 1> current_multi_index = {1};
  multi_index_handler_1D_->SetCurrentIndex(current_multi_index);
  ASSERT_THAT(multi_index_handler_1D_->CollapseDimension(Dimension{0}), Eq(0));

}

class MultiIndexHandler2D : public Test {
 public:
  MultiIndexHandler2D() {
    std::array<int, 2> last_knot_offset_2D = {10, 3};
    multi_index_handler_2D_ = std::make_unique<util::MultiIndexHandler<2>>(last_knot_offset_2D);
  }

 protected:
  std::unique_ptr<util::MultiIndexHandler<2>> multi_index_handler_2D_;
};

TEST_F(MultiIndexHandler2D, Returns1DIndex0AfterConstruction) { // NOLINT
  ASSERT_THAT(multi_index_handler_2D_->GetCurrent1DIndex(), Eq(0));
}

TEST_F(MultiIndexHandler2D, Returns1DIndex23AfterSettingCurrentIndexTo3And2) { // NOLINT
  std::array<int, 2> current_2D_index = {3, 2};
  multi_index_handler_2D_->SetCurrentIndex(current_2D_index);
  ASSERT_THAT(multi_index_handler_2D_->GetCurrent1DIndex(), Eq(23));
}

TEST_F(MultiIndexHandler2D, Returns2DIndex3And2AfterSettingCurrent1DIndexTo23) { // NOLINT
  multi_index_handler_2D_->SetCurrentIndex(23);
  ASSERT_THAT(multi_index_handler_2D_->GetCurrentIndex()[0], Eq(3));
  ASSERT_THAT(multi_index_handler_2D_->GetCurrentIndex()[1], Eq(2));
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
