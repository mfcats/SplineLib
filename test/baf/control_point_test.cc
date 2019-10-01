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

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

using namespace splinelib::src;

class AControlPoint : public Test {
 public:
  AControlPoint() : control_point({1.0, 2.0}), control_point_b({2.0, -1.0}), control_point_c({1.0, 0.7, 1.9}),
      transMatrix({std::array<double, 4>({0.5, 0.0, 0.866, 0.0}),
                   std::array<double, 4>({0.0, 1.0, 0.0, 0.0}),
                   std::array<double, 4>({-0.866, 0.0, 0.5, 0.0}),
                   std::array<double, 4>({0.0, 0.0, 0.0, 1.0})}),
      scaling({1.0, 1.0, 1.0}) {}

 protected:
  baf::ControlPoint control_point;
  baf::ControlPoint control_point_b;
  baf::ControlPoint control_point_c;
  std::array<std::array<double, 4>, 4> transMatrix;
  std::array<double, 3> scaling;
};

TEST_F(AControlPoint, ReturnsCorrectDimension) { // NOLINT
  ASSERT_THAT(control_point.GetDimension(), 2);
}

TEST_F(AControlPoint, Returns1For0Dimension) { // NOLINT
  ASSERT_THAT(control_point.GetValue(0), DoubleEq(1.0));
}

TEST_F(AControlPoint, Returns2For1Dimension) { // NOLINT
  ASSERT_THAT(control_point.GetValue(1), DoubleEq(2.0));
}

TEST_F(AControlPoint, Returns3ForSum0Dimension) { // NOLINT
  ASSERT_THAT((control_point + control_point_b).GetValue(0), DoubleEq(3.0));
}

TEST_F(AControlPoint, PerformsTransformationCorrectly) { // NOLINT
  baf::ControlPoint control_point_t = control_point_c.Transform(transMatrix, scaling);
  ASSERT_THAT(control_point_t.GetValue(0), DoubleEq(2.1454));
}

TEST_F(AControlPoint, PerformsScalingTransformation) { // NOLINT
  std::array<double, 3> customScaling = {0.5, 1, 0.7};
  baf::ControlPoint control_point_t = control_point_c.Transform(transMatrix, customScaling);
  ASSERT_THAT(control_point_t.GetValue(0), DoubleNear(1.40178, 0.00001));
}
