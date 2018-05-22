/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "gmock/gmock.h"

#include "multi_index_handler.h"

using testing::Test;
using testing::Eq;

class MultiHandlers : public Test {
 public:
  MultiHandlers() {
    std::array<int, 1> lastKnotOffset1D = {10};
    //currentIndex[0] = 10;
    multiIndexHandler1D = std::make_unique<util::MultiIndexHandler<1>>(lastKnotOffset1D);
    std::array<int, 2> lastKnotOffset2D = {10, 3};
    //lastKnotOffset2D[0] = 10;
    //lastKnotOffset2D[1] = 3;
    multiIndexHandler2D = std::make_unique<util::MultiIndexHandler<2>>(lastKnotOffset2D);
    std::array<int, 3> lastKnotOffset3D = {4, 3, 5};
    //lastKnotOffset3D[0] = 4;
    //lastKnotOffset3D[1] = 3;
    //lastKnotOffset3D[2] = 5;
    multiIndexHandler3D = std::make_unique<util::MultiIndexHandler<3>>(lastKnotOffset3D);
  }
 protected:
  std::unique_ptr<util::MultiIndexHandler<1>> multiIndexHandler1D;
  std::unique_ptr<util::MultiIndexHandler<2>> multiIndexHandler2D;
  std::unique_ptr<util::MultiIndexHandler<3>> multiIndexHandler3D;
};

TEST_F(MultiHandlers, 1D_multiIndex_ind0) {
  ASSERT_THAT(multiIndexHandler1D->Get1DIndex(), Eq(0));
}

TEST_F(MultiHandlers, 1D_multiIndex_ind5) {
  std::array<int, 1> currentIndex = {5};
  //currentIndex[0] = 5;
  multiIndexHandler1D->SetIndices(currentIndex);
  ASSERT_THAT(multiIndexHandler1D->Get1DIndex(), Eq(5));
}

TEST_F(MultiHandlers, 2D_multiIndex_ind0) {
  ASSERT_THAT(multiIndexHandler2D->Get1DIndex(), Eq(0));
}

TEST_F(MultiHandlers, 2D_multiIndex_ind3_2) {
  std::array<int, 2> currentIndex = {3, 2};
  //currentIndex[0] = 3;
  //currentIndex[1] = 2;
  multiIndexHandler2D->SetIndices(currentIndex);
  ASSERT_THAT(multiIndexHandler2D->Get1DIndex(), Eq(23));
}

TEST_F(MultiHandlers, 3D_multiIndex_ind0) {
  ASSERT_THAT(multiIndexHandler3D->Get1DIndex(), Eq(0));
}

TEST_F(MultiHandlers, 3D_multiIndex_ind2_1_2) {
  std::array<int, 3> currentIndex = {2, 1, 2};
  //currentIndex[0] = 2;
  //currentIndex[1] = 1;
  //currentIndex[2] = 2;
  multiIndexHandler3D->SetIndices(currentIndex);
  ASSERT_THAT(multiIndexHandler3D->Get1DIndex(), Eq(30));
}
