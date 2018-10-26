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

#include "b_spline.h"
#include "nurbs.h"
#include "connectivity_handler.h"

using testing::Test;
using testing::DoubleEq;

class IGA2D : public Test {
 public:
  std::array<baf::KnotVector, 2> knot_vector =
      {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.4}, ParamCoord{0.5},
                        ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}),
       baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{0.5},
                        ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})};

  std::array<Degree, 2> degree = {Degree{3}, Degree{3}};

  std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1};

  std::vector<baf::ControlPoint> control_points = {
      baf::ControlPoint(std::vector<double>({0.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.16, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.32, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.48, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.64, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.8, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -1.0, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.16, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.32, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.48, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.64, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.8, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -0.66, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.16, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.32, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.48, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.64, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.8, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -0.33, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.16, 0.16, 0.0})),
      baf::ControlPoint(std::vector<double>({0.32, 0.32, 0.0})),
      baf::ControlPoint(std::vector<double>({0.48, 0.48, 0.0})),
      baf::ControlPoint(std::vector<double>({0.64, 0.64, 0.0})),
      baf::ControlPoint(std::vector<double>({0.8, 0.8, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0})),

      baf::ControlPoint(std::vector<double>({-0.33, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.16, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.32, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.48, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.64, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.8, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 1.0, 0.0})),

      baf::ControlPoint(std::vector<double>({-0.66, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.16, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.32, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.48, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.64, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.8, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 1.0, 0.0})),

      baf::ControlPoint(std::vector<double>({-1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.16, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.32, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.48, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.64, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.8, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 1.0, 0.0}))
  };

  std::array<std::shared_ptr<baf::KnotVector>, 2> kv_ptr = {std::make_shared<baf::KnotVector>(knot_vector[0]),
                                                            std::make_shared<baf::KnotVector>(knot_vector[1])};
  std::shared_ptr<spl::NURBS<2>> nurbs_ = std::make_shared<spl::NURBS<2>>(kv_ptr, degree, control_points, weights);
};

TEST_F(IGA2D, TestConnectivityHandler) { // NOLINT
  iga::ConnectivityHandler connectivity_handler = iga::ConnectivityHandler(nurbs_);
  std::vector<int> e1 = {1, 2, 3, 4, 8, 9, 10, 11, 15, 16, 17, 18, 22, 23, 24, 25};
  std::vector<int> e2 = {2, 3, 4, 5, 9, 10, 11, 12, 16, 17, 18, 19, 23, 24, 25, 26};
  std::vector<int> e3 = {3, 4, 5, 6, 10, 11, 12, 13, 17, 18, 19, 20, 24, 25, 26, 27};
  std::vector<int> e4 = {4, 5, 6, 7, 11, 12, 13, 14, 18, 19, 20, 21, 25, 26, 27, 28};
  std::vector<int> e5 = {22, 23, 24, 25, 29, 30, 31, 32, 36, 37, 38, 39, 43, 44, 45, 46};
  std::vector<int> e6 = {23, 24, 25, 26, 30, 31, 32, 33, 37, 38, 39, 40, 44, 45, 46, 47};
  std::vector<int> e7 = {24, 25, 26, 27, 31, 32, 33, 34, 38, 39, 40, 41, 45, 46, 47, 48};
  std::vector<int> e8 = {25, 26, 27, 28, 32, 33, 34, 35, 39, 40, 41, 42, 46, 47, 48, 49};

  std::vector<std::vector<int>> connectivity_matlab = {e1, e2, e3, e4, e5, e6, e7, e8};
  std::vector<std::vector<int>> connectivity_splinelib = connectivity_handler.GetConnectivity();
  ASSERT_EQ(connectivity_matlab.size(), connectivity_splinelib.size());
  for (uint64_t i = 0; i < connectivity_matlab.size(); ++i) {
    ASSERT_EQ(connectivity_matlab[i].size(), connectivity_splinelib[i].size());
  }
  for (uint64_t i = 0; i < connectivity_matlab.size(); ++i) {
    for (uint64_t j = 0; j < connectivity_matlab[i].size(); ++j) {
      ASSERT_EQ(connectivity_matlab[i][j], connectivity_splinelib[i][j]);
    }
  }
}
