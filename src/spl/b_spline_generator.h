/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_B_SPLINE_GENERATOR_H_
#define SRC_SPL_B_SPLINE_GENERATOR_H_

#include <vector>

#include "physical_space.h"
#include "spline_generator.h"

namespace spl {
template<int DIM>
class BSplineGenerator : public SplineGenerator<DIM> {
 public:
  BSplineGenerator() = default;
  virtual ~BSplineGenerator() = default;

  BSplineGenerator(std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vector,
                   std::array<int, DIM> degree,
                   const std::vector<baf::ControlPoint> &control_points) {
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = knot_vector[i]->GetNumberOfKnots() - degree[i] - 1;
    }

    std::array<baf::KnotVector, DIM> knot_vectors;
    for (int dim = 0; dim < DIM; ++dim) {
      knot_vectors[dim] = *(knot_vector[dim]);
    }

    this->physical_space_ = std::make_shared<PhysicalSpace<DIM>>(control_points, number_of_points);
    this->parameter_space_ = std::make_shared<ParameterSpace<DIM>>(knot_vectors, degree);
  }

  BSplineGenerator(PhysicalSpace<DIM> physical_space, ParameterSpace<DIM> parameter_space) {
    this->physical_space_ = std::make_shared<PhysicalSpace<DIM>>(physical_space);
    this->parameter_space_ = std::make_shared<ParameterSpace<DIM>>(parameter_space);
  }

  std::shared_ptr<PhysicalSpace<DIM>> GetPhysicalSpace() {
    return physical_space_;
  }

 protected:
  std::shared_ptr<PhysicalSpace<DIM>> physical_space_;
};
}  // namespace spl

#endif  // SRC_SPL_B_SPLINE_GENERATOR_H_