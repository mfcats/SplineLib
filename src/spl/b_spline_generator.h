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

namespace splinelib::src::spl {
template<int PARAMETRIC_DIMENSIONALITY>
class BSplineGenerator : public SplineGenerator<PARAMETRIC_DIMENSIONALITY> {
 public:
  BSplineGenerator() = default;
  BSplineGenerator(const BSplineGenerator<PARAMETRIC_DIMENSIONALITY> &other) = default;
  BSplineGenerator(BSplineGenerator<PARAMETRIC_DIMENSIONALITY> &&other) = delete;
  BSplineGenerator & operator=(const BSplineGenerator<PARAMETRIC_DIMENSIONALITY> &rhs) = default;
  BSplineGenerator & operator=(BSplineGenerator<PARAMETRIC_DIMENSIONALITY> &&rhs) = delete;
  ~BSplineGenerator() override = default;

  BSplineGenerator(baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> knot_vector,
      std::array<Degree, PARAMETRIC_DIMENSIONALITY> degree, const std::vector<baf::ControlPoint> &control_points)
      : SplineGenerator<PARAMETRIC_DIMENSIONALITY>(knot_vector, degree) {
    std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points;
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      number_of_points[i] = knot_vector[i]->GetNumberOfKnots() - degree[i].Get() - 1;
    }
    physical_space_ = std::make_shared<PhysicalSpace<PARAMETRIC_DIMENSIONALITY>>(control_points, number_of_points);
  }

  BSplineGenerator(std::shared_ptr<PhysicalSpace<PARAMETRIC_DIMENSIONALITY>> physical_space,
                   std::shared_ptr<ParameterSpace<PARAMETRIC_DIMENSIONALITY>> parameter_space) {
    this->parameter_space_ = parameter_space;
    physical_space_ = physical_space;
  }

  std::shared_ptr<PhysicalSpace<PARAMETRIC_DIMENSIONALITY>> GetPhysicalSpace() const {
    return physical_space_;
  }

 protected:
  std::shared_ptr<PhysicalSpace<PARAMETRIC_DIMENSIONALITY>> physical_space_;
};
}  //  namespace splinelib::src::spl

#endif  //  SRC_SPL_B_SPLINE_GENERATOR_H_
