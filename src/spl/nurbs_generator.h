/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_NURBS_GENERATOR_H_
#define SRC_SPL_NURBS_GENERATOR_H_

#include <vector>

#include "spline_generator.h"
#include "weighted_physical_space.h"

namespace spl {
template<int DIM>
class NURBSGenerator : public SplineGenerator<DIM> {
 public:
  NURBSGenerator() = default;
  virtual ~NURBSGenerator() = default;

  NURBSGenerator(std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vector,
                 std::array<int, DIM> degree,
                 const std::vector<baf::ControlPoint> &control_points, std::vector<double> weights) {
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = knot_vector[i]->GetNumberOfKnots() - degree[i] - 1;
    }

    std::array<baf::KnotVector, DIM> knot_vectors;
    for (int dim = 0; dim < DIM; ++dim) {
      knot_vectors[dim] = *(knot_vector[dim]);
    }


    this->physical_space_ptr = std::make_shared<WeightedPhysicalSpace<DIM>>(control_points, weights, number_of_points);
    this->parameter_space_ptr = std::make_shared<ParameterSpace<DIM>>(knot_vectors, degree);
  }

  NURBSGenerator(WeightedPhysicalSpace<DIM> weigthed_physical_space, ParameterSpace<DIM> parameter_space) {
    this->physical_space_ptr = std::make_shared<WeightedPhysicalSpace<DIM>>(weigthed_physical_space);
    this->parameter_space_ptr = std::make_shared<ParameterSpace<DIM>>(parameter_space);
  }

  std::shared_ptr<WeightedPhysicalSpace<DIM>> GetWeightedPhysicalSpace() {
    return physical_space_ptr;
  }

 protected:
  std::shared_ptr<WeightedPhysicalSpace<DIM>> physical_space_ptr;
};
}  // namespace spl

#endif  // SRC_SPL_NURBS_GENERATOR_H_
