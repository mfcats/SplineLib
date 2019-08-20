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

#include "control_point.h"
#include "bezier_segment.h"

using testing::Test;

class BezierSegmentTest : public Test {  // NOLINT
 public:
  BezierSegmentTest() {
    std::array<Degree, 1> degree = {Degree{4}};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({3.0, 1.0})),
        baf::ControlPoint(std::vector<double>({4.0, 0.0}))
    };
    std::array<bool, 1> is_bezier_in_direction = {true};
    std::array<int, 1> num_control_points = {5};
    bezier_segment_ = std::make_shared<spl::BezierSegment<1>>(degree, control_points, is_bezier_in_direction,
        num_control_points);
  }

 protected:
  std::shared_ptr<spl::BezierSegment<1>> bezier_segment_;
};

TEST_F(BezierSegmentTest, ReducesDegreeFrom5To4Correctly) {  // NOLINT

  // TODO: Create proper test cases!
  bezier_segment_->ReduceDegree(0);
}
