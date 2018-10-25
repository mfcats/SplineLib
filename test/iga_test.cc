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
      {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.25}, ParamCoord{0.5}, ParamCoord{0.5},
                        ParamCoord{0.75}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}),
       baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.25}, ParamCoord{0.5}, ParamCoord{0.5},
                        ParamCoord{0.75}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})};
  std::array<Degree, 2> degree = {Degree{2}, Degree{2}};
  std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<baf::ControlPoint> control_points = {
      baf::ControlPoint(std::vector<double>({-1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 1.0, 0.0}))
  };
  std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector_ptr =
      {std::make_shared<baf::KnotVector>(knot_vector[0]),
       std::make_shared<baf::KnotVector>(knot_vector[1])};
  std::shared_ptr<spl::NURBS<2>> nurbs_ = std::make_shared<spl::NURBS<2>>(knot_vector_ptr, degree, control_points,
      weights);
};

TEST_F(IGA2D, Test1) { // NOLINT
  iga::ConnectivityHandler connectivity_handler = iga::ConnectivityHandler(nurbs_);

}

