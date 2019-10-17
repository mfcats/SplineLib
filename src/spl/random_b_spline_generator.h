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

#include <memory>
#include <vector>

#include "b_spline_generator.h"
#include "random_spline_utils.h"

namespace splinelib::src::spl {
template<int PARAMETRIC_DIMENSIONALITY>
class RandomBSplineGenerator : public BSplineGenerator<PARAMETRIC_DIMENSIONALITY> {
 public:
  RandomBSplineGenerator() = default;
  RandomBSplineGenerator(const RandomBSplineGenerator<PARAMETRIC_DIMENSIONALITY> &other) = delete;
  RandomBSplineGenerator(RandomBSplineGenerator<PARAMETRIC_DIMENSIONALITY> &&other) = delete;
  RandomBSplineGenerator & operator=(const RandomBSplineGenerator<PARAMETRIC_DIMENSIONALITY> &rhs) = delete;
  RandomBSplineGenerator & operator=(RandomBSplineGenerator<PARAMETRIC_DIMENSIONALITY> &&rhs) = delete;
  ~RandomBSplineGenerator() override = default;

  RandomBSplineGenerator(std::array<ParametricCoordinate, 2> coord_limits, int max_degree, int dimension) {
    std::array<Degree, PARAMETRIC_DIMENSIONALITY> degrees =
        RandomSplineUtils<PARAMETRIC_DIMENSIONALITY>::GetRandomDegrees(max_degree);
    baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> knot_vectors =
        RandomSplineUtils<PARAMETRIC_DIMENSIONALITY>::GetRandomKnotVectors(coord_limits, degrees);
    std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points =
        RandomSplineUtils<PARAMETRIC_DIMENSIONALITY>::GetNumberOfPoints(degrees, knot_vectors);
    std::vector<baf::ControlPoint> control_points =
        RandomSplineUtils<PARAMETRIC_DIMENSIONALITY>::GetRandomControlPoints(dimension, number_of_points);
    this->physical_space_ =
        std::make_shared<PhysicalSpace<PARAMETRIC_DIMENSIONALITY>>(control_points, number_of_points);
    this->parameter_space_ = std::make_shared<ParameterSpace<PARAMETRIC_DIMENSIONALITY>>(knot_vectors, degrees);
  }
};
}  // namespace splinelib::src::spl

#endif  // SRC_SPL_RANDOM_B_SPLINE_GENERATOR_H_
