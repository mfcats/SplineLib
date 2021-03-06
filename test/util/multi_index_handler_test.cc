/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "multi_index_handler.h"

#include "gmock/gmock.h"

using testing::Test;
using testing::Eq;

class MultiHandler1D : public Test {
 public:
  MultiHandler1D() {
    std::array<int, 1> lastKnotOffset1D = {10};
    multiIndexHandler1D = std::make_unique<util::MultiIndexHandler<1>>(lastKnotOffset1D);
  }

 protected:
  std::unique_ptr<util::MultiIndexHandler<1>> multiIndexHandler1D;
};

TEST_F(MultiHandler1D, Returns1DIndex0AfterConstruction) { // NOLINT
  ASSERT_THAT(multiIndexHandler1D->Get1DIndex(), Eq(0));
}

TEST_F(MultiHandler1D, Returns1DIndex5AfterSettingCurrentIndexTo5) { // NOLINT
  std::array<int, 1> currentIndex = {5};
  multiIndexHandler1D->SetIndices(currentIndex);
  ASSERT_THAT(multiIndexHandler1D->Get1DIndex(), Eq(5));
}

class MultiHandler2D : public Test {
 public:
  MultiHandler2D() {
    std::array<int, 2> lastKnotOffset2D = {10, 3};
    multiIndexHandler2D = std::make_unique<util::MultiIndexHandler<2>>(lastKnotOffset2D);
  }

 protected:
  std::unique_ptr<util::MultiIndexHandler<2>> multiIndexHandler2D;
};

TEST_F(MultiHandler2D, Returns1DIndex0AfterConstruction) { // NOLINT
  ASSERT_THAT(multiIndexHandler2D->Get1DIndex(), Eq(0));
}

TEST_F(MultiHandler2D, Returns1DIndex23AfterSettingCurrentIndexTo3And2) { // NOLINT
  std::array<int, 2> currentIndex = {3, 2};
  multiIndexHandler2D->SetIndices(currentIndex);
  ASSERT_THAT(multiIndexHandler2D->Get1DIndex(), Eq(23));
}

TEST_F(MultiHandler2D, Returns2DIndex3And2AfterSettingCurrent1DIndexTo23) { // NOLINT
  multiIndexHandler2D->Set1DIndex(23);
  ASSERT_THAT(multiIndexHandler2D->GetIndices()[0], Eq(3));
  ASSERT_THAT(multiIndexHandler2D->GetIndices()[1], Eq(2));
}

class MultiHandler3D : public Test {
 public:
  MultiHandler3D() {
    std::array<int, 3> lastKnotOffset3D = {4, 3, 5};
    multiIndexHandler3D = std::make_unique<util::MultiIndexHandler<3>>(lastKnotOffset3D);
  }
 protected:
  std::unique_ptr<util::MultiIndexHandler<3>> multiIndexHandler3D;
};

TEST_F(MultiHandler3D, Returns1DIndex0AfterConstruction) { // NOLINT
  ASSERT_THAT(multiIndexHandler3D->Get1DIndex(), Eq(0));
}

TEST_F(MultiHandler3D, Returns1DIndex30AfterSettingCurrentIndexTo2And1And2) { // NOLINT
  std::array<int, 3> currentIndex = {2, 1, 2};
  multiIndexHandler3D->SetIndices(currentIndex);
  ASSERT_THAT(multiIndexHandler3D->Get1DIndex(), Eq(30));
}

TEST_F(MultiHandler3D, ReturnsIndex3And1And0AfterAdding1DIndex7) { // NOLINT
  *multiIndexHandler3D + 7;
  ASSERT_THAT(multiIndexHandler3D->Get1DIndex(), Eq(7));
  ASSERT_THAT(multiIndexHandler3D->GetIndices()[0], Eq(3));
  ASSERT_THAT(multiIndexHandler3D->GetIndices()[1], Eq(1));
  ASSERT_THAT(multiIndexHandler3D->GetIndices()[2], Eq(0));
}

TEST_F(MultiHandler3D, ReturnsDifferenceIndex0And1And4AfterAdding1DIndex7) { // NOLINT
  *multiIndexHandler3D + 7;
  ASSERT_THAT(multiIndexHandler3D->Get1DIndex(), Eq(7));
  ASSERT_THAT(multiIndexHandler3D->GetDifferenceIndices()[0], Eq(0));
  ASSERT_THAT(multiIndexHandler3D->GetDifferenceIndices()[1], Eq(1));
  ASSERT_THAT(multiIndexHandler3D->GetDifferenceIndices()[2], Eq(4));
}

TEST_F(MultiHandler3D, ReturnsCorrectLengthOf60) { // NOLINT
  ASSERT_THAT(multiIndexHandler3D->Get1DLength(), Eq(60));
}

TEST_F(MultiHandler3D, ReturnsCorrectLengthAfterExtractingDimension1) { // NOLINT
  std::array<int, 3> indices = {2, 1, 3};
  multiIndexHandler3D->SetIndices(indices);
  ASSERT_THAT(multiIndexHandler3D->ExtractDimension(1), Eq(14));
}
