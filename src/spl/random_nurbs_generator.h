/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_RANDOM_NURBS_GENERATOR_H_
#define SRC_SPL_RANDOM_NURBS_GENERATOR_H_

#include <memory>
#include <vector>

#include "nurbs_generator.h"
#include "random_spline_utils.h"

namespace spl {
template<int DIM>
class RandomNURBSGenerator : public NURBSGenerator<DIM> {
 public:
  RandomNURBSGenerator() = default;
  virtual ~RandomNURBSGenerator() = default;

  RandomNURBSGenerator(std::array<ParamCoord, 2> coord_limits, int max_degree, int dimension) {
    std::array<Degree, DIM> degrees = RandomSplineUtils<DIM>::GetRandomDegrees(max_degree);
    KnotVectors<DIM> knot_vectors = RandomSplineUtils<DIM>::GetRandomKnotVectors(coord_limits, degrees);
    std::array<int, DIM> number_of_points = RandomSplineUtils<DIM>::GetNumberOfPoints(degrees, knot_vectors);
    std::vector<baf::ControlPoint>
        control_points = RandomSplineUtils<DIM>::GetRandomControlPoints(dimension, number_of_points);
    std::vector<double> weights = RandomSplineUtils<DIM>::GetRandomWeights(number_of_points);
    this->physical_space_ = std::make_shared<WeightedPhysicalSpace<DIM>>(control_points, weights, number_of_points);
    this->parameter_space_ = std::make_shared<ParameterSpace<DIM>>(knot_vectors, degrees);
  }
};
}  // namespace spl

#endif  // SRC_SPL_RANDOM_NURBS_GENERATOR_H_
