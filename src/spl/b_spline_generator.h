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

#include <memory>
#include <vector>

#include "physical_space.h"
#include "spline_generator.h"

namespace spl {
template<int DIM>
class BSplineGenerator : public SplineGenerator<DIM> {
 public:
  BSplineGenerator() = default;
  virtual ~BSplineGenerator() = default;

  BSplineGenerator(KnotVectors<DIM> knot_vector, std::array<Degree, DIM> degree,
                   const std::vector<baf::ControlPoint> &control_points) : SplineGenerator<DIM>(knot_vector, degree) {
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = knot_vector[i]->GetNumberOfKnots() - degree[i].get() - 1;
    }
    physical_space_ = std::make_shared<PhysicalSpace<DIM>>(control_points, number_of_points);
  }

  BSplineGenerator(std::shared_ptr<PhysicalSpace<DIM>> physical_space,
                   std::shared_ptr<ParameterSpace<DIM>> parameter_space) {
    this->parameter_space_ = parameter_space;
    physical_space_ = physical_space;
  }

  std::shared_ptr<PhysicalSpace<DIM>> GetPhysicalSpace() const {
    return physical_space_;
  }

 protected:
  std::shared_ptr<PhysicalSpace<DIM>> physical_space_;
};
}  //  namespace spl

#endif  //  SRC_SPL_B_SPLINE_GENERATOR_H_
