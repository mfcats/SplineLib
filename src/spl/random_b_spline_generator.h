/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_RANDOM_B_SPLINE_GENERATOR_H_
#define SRC_SPL_RANDOM_B_SPLINE_GENERATOR_H_

#include <iostream>
#include <memory>
#include <vector>

#include "physical_space.h"
#include "b_spline_generator.h"

namespace spl {
template<int DIM>
class RandomBSplineGenerator : public BSplineGenerator<DIM> {
 public:
  RandomBSplineGenerator() = default;
  virtual ~RandomBSplineGenerator() = default;

  RandomBSplineGenerator(std::array<std::array<ParamCoord, 2>, DIM> param_coord_limits,
                         std::array<Degree, DIM> degree,
                         int dimension) {
    std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vector;
    std::array<std::vector<ParamCoord>, DIM> param_coord_vectors;
    for (int i = 0; i < DIM; ++i) {
      for (int j = 0; j < degree[i].get() + 1; ++j) {
        param_coord_vectors[i].push_back(param_coord_limits[i][0]);
      }
      for (int j = 0; j < degree[i].get() + 1; ++j) {
        param_coord_vectors[i].push_back(param_coord_limits[i][1]);
      }
      baf::KnotVector knt_vec(param_coord_vectors[i]);
      knot_vector[i] = std::make_shared<baf::KnotVector>(knt_vec);
    }
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = degree[i].get() + 1;
    }
    std::vector<baf::ControlPoint> control_points;
    util::MultiIndexHandler<DIM> point_handler(number_of_points);
    for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
      std::vector<double> coordinates;
      for (int j = 0; j < dimension; ++j) {
        coordinates.push_back((rand() % 100) / 10.0 - 5.0);
      }
      control_points.emplace_back(coordinates);
    }
    this->physical_space_ = std::make_shared<PhysicalSpace<DIM>>(control_points, number_of_points);
    this->parameter_space_ = std::make_shared<ParameterSpace<DIM>>(knot_vector, degree);
  }
};
}  //  namespace spl

#endif  //  SRC_SPL_RANDOM_B_SPLINE_GENERATOR_H_
