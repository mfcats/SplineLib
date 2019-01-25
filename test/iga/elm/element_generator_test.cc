/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <b_spline.h>
#include "gmock/gmock.h"

#include "element_generator.h"

using testing::Test;

class A1DElementGenerator : public Test {
 public:
  std::array<baf::KnotVector, 1> knot_vector =
      {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                        ParamCoord{2}, ParamCoord{3}, ParamCoord{4}, ParamCoord{4},
                        ParamCoord{5}, ParamCoord{5}, ParamCoord{5}})};
  std::array<Degree, 1> degree = {Degree{2}};
  std::vector<baf::ControlPoint> control_points = {
      baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({2.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({3.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({4.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({5.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({6.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({7.0, 0.0, 0.0}))};
  KnotVectors<1> kv_ptr = {std::make_shared<baf::KnotVector>(knot_vector[0])};
  std::shared_ptr<spl::BSpline<1>> b_spline = std::make_shared<spl::BSpline<1>>(kv_ptr, degree, control_points);

  A1DElementGenerator() : element_generator(b_spline) {}

 protected:
  iga::elm::ElementGenerator<1> element_generator;
};

TEST_F(A1DElementGenerator, ReturnsCorrectNumberOfElements) { // NOLINT
  ASSERT_THAT(element_generator.GetElementList(0).size(), 5);
}

TEST_F(A1DElementGenerator, ReturnsElementsWithCorrectNodes) { // NOLINT
  auto element_list = element_generator.GetElementList(0);
  for (auto element = 0u; element < element_list.size(); element++) {
    ASSERT_THAT(element_list[element].GetLowerBound().get(), element);
    ASSERT_THAT(element_list[element].GetUpperBound().get(), element + 1);
  }
}
