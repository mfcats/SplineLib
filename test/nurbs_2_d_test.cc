/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <numeric_settings.h>
#include "gmock/gmock.h"

#include "nurbs.h"

using testing::Test;
using testing::DoubleNear;

class NurbsEx4_3 : public Test {
 public:
  NurbsEx4_3() {

    std::array<baf::KnotVector, 2> knot_vector = {baf::KnotVector({0, 0, 0, 1, 1, 1}),
                                                  baf::KnotVector({0, 0, 0, 1, 1, 1})};
    std::array<int, 2> degree = {2, 2};
    std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1, 2, 1};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({3.0, 0.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 0.5})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({2.5, 3.5})),
        baf::ControlPoint(std::vector<double>({5.0, 2.0}))
    };
    nurbs_ = std::make_unique<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
  }

 protected:
  std::unique_ptr<spl::NURBS<2>> nurbs_;
};

TEST_F(NurbsEx4_3, Returns2_0For2_5And1_0AndDim0) {
  ASSERT_THAT(nurbs_->Evaluate({0.5, 1.0}, {0})[0], DoubleNear(2.5, util::NumericSettings<double>::kEpsilon()));
}

TEST_F(NurbsEx4_3, Returns3_63For2_5And1_0AndDim1) {
  ASSERT_THAT(nurbs_->Evaluate({0.5, 1.0}, {1})[0], DoubleNear(3.0, util::NumericSettings<double>::kEpsilon()));
}
