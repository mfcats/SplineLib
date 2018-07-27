/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_B_SPLINE_GENERATOR_H
#define SPLINELIB_B_SPLINE_GENERATOR_H

#include "spline_generator.h"
#include <memory>

namespace spl {
template<int DIM>
class BSplineGenerator : public SplineGenerator<DIM> {
 public:
  BSplineGenerator(std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vector,
                   std::array<int, DIM> degree,
                   const std::vector<baf::ControlPoint> &control_points) : SplineGenerator<DIM>(knot_vector, degree) {
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = knot_vector[i]->GetNumberOfKnots() - degree[i] - 1;
    }
    physical_space_ = PhysicalSpace<DIM>(control_points, number_of_points);
  }

  BSplineGenerator(std::shared_ptr<PhysicalSpace<1>> physical_space, std::shared_ptr<ParameterSpace<1>> parameter_space) {
    this->parameter_space_ = parameter_space;
    physical_space_ = physical_space;
  }

  std::shared_ptr<PhysicalSpace<DIM>> GetPhysicalSpace() const {
    return physical_space_;
  }

 private:
  std::shared_ptr<PhysicalSpace<DIM>> physical_space_;
};
}

#endif //SPLINELIB_B_SPLINE_GENERATOR_H
