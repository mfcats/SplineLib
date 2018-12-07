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
#include "connectivity_handler.h"
#include "test_spline.h"

using testing::Eq;

TEST_F(AnIGATestSpline, TestConnectivityHandler) { // NOLINT
  iga::ConnectivityHandler<2> connectivity_handler = iga::ConnectivityHandler<2>(nurbs_);
  std::vector<int> e1 = {1, 2, 3, 4, 8, 9, 10, 11, 15, 16, 17, 18, 22, 23, 24, 25};
  std::vector<int> e2 = {2, 3, 4, 5, 9, 10, 11, 12, 16, 17, 18, 19, 23, 24, 25, 26};
  std::vector<int> e3 = {3, 4, 5, 6, 10, 11, 12, 13, 17, 18, 19, 20, 24, 25, 26, 27};
  std::vector<int> e4 = {4, 5, 6, 7, 11, 12, 13, 14, 18, 19, 20, 21, 25, 26, 27, 28};
  std::vector<int> e5 = {22, 23, 24, 25, 29, 30, 31, 32, 36, 37, 38, 39, 43, 44, 45, 46};
  std::vector<int> e6 = {23, 24, 25, 26, 30, 31, 32, 33, 37, 38, 39, 40, 44, 45, 46, 47};
  std::vector<int> e7 = {24, 25, 26, 27, 31, 32, 33, 34, 38, 39, 40, 41, 45, 46, 47, 48};
  std::vector<int> e8 = {25, 26, 27, 28, 32, 33, 34, 35, 39, 40, 41, 42, 46, 47, 48, 49};
  std::vector<std::vector<int>> connectivity_matlab = {e1, e2, e3, e4, e5, e6, e7, e8};
  for (uint64_t i = 0; i < connectivity_matlab.size(); ++i) {
    for (uint64_t j = 0; j < connectivity_matlab[i].size(); ++j) {
      ASSERT_THAT(connectivity_handler.GetGlobalIndex(i, j), Eq(connectivity_matlab[i][j]));
    }
  }
}

TEST_F(AnIGATestSpline3, TestConnectivityHandler) { // NOLINT
  iga::ConnectivityHandler<3> connectivity_handler = iga::ConnectivityHandler<3>(nurbs_);
     ASSERT_THAT(connectivity_handler.GetGlobalIndex(15, 10), Eq(188));
}