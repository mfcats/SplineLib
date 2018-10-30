/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "surface_generator.h"

spl::SurfaceGenerator::SurfaceGenerator(std::shared_ptr<spl::NURBS<1>> const nurbs_T,
                                        std::shared_ptr<spl::NURBS<1>> const nurbs_C) {
  this->parameter_space_= JoinParameterSpaces(nurbs_T->GetParameterSpace(), nurbs_C->GetParameterSpace());
  this->physical_space_ = JoinPhysicalSpaces(nurbs_T->GetPhysicalSpace(), nurbs_C->GetPhysicalSpace());
}

std::shared_ptr<spl::ParameterSpace<2>> spl::SurfaceGenerator::JoinParameterSpaces(
    std::shared_ptr<ParameterSpace<1>> const space_1, std::shared_ptr<ParameterSpace<1>> const space_2) const {
  std::array<std::shared_ptr<baf::KnotVector>, 2> joined_knot_vector =
      {space_1->GetKnotVector(0), space_2->GetKnotVector(0)};
  std::array<Degree, 2> joined_degree = {space_1->GetDegree(0), space_2->GetDegree(0)};
  return std::make_shared<ParameterSpace<2>>(ParameterSpace<2>(joined_knot_vector, joined_degree));
}

std::shared_ptr<spl::WeightedPhysicalSpace<2>> spl::SurfaceGenerator::JoinPhysicalSpaces(
    std::shared_ptr<spl::PhysicalSpace<1>> const space_1, std::shared_ptr<spl::PhysicalSpace<1>> const space_2) const {
  std::array<int, 2> joined_number_of_points =
      {space_1->GetNumberOfControlPoints(), space_2->GetNumberOfControlPoints()};
  std::vector<baf::ControlPoint> joined_control_points;
  std::vector<double> joined_weights;
  for (int i = 0; i < space_2->GetNumberOfControlPoints(); ++i) {
    std::array<int, 1> index_space_2 = {i};
    for (int j = 0; j < space_1->GetNumberOfControlPoints(); ++j) {
      std::array<int, 1> index_space_1 = {j};
      joined_control_points.emplace_back(space_1->GetControlPoint(index_space_1) +
                                         space_2->GetControlPoint(index_space_2));
      joined_weights.emplace_back(space_1->GetWeight(index_space_1) * space_2->GetWeight(index_space_2));
    }
  }
  return std::make_shared<spl::WeightedPhysicalSpace<2>>(spl::WeightedPhysicalSpace<2>(
      joined_control_points, joined_weights, joined_number_of_points));
}


