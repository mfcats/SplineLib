/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef TEST_SPL_RANDOM_RANDOM_NURBS_GENERATOR_H_
#define TEST_SPL_RANDOM_RANDOM_NURBS_GENERATOR_H_

#include <memory>
#include <vector>

#include "src/spl/nurbs_generator.h"
#include "random_spline_utils.h"

namespace splinelib::test::spl::random {
using namespace splinelib::src;
template<int PARAMETRIC_DIMENSIONALITY>
 class RandomNURBSGenerator : public src::spl::NURBSGenerator<PARAMETRIC_DIMENSIONALITY> {
 public:
  RandomNURBSGenerator() = default;
  RandomNURBSGenerator(const RandomNURBSGenerator<PARAMETRIC_DIMENSIONALITY> &other) = delete;
  RandomNURBSGenerator(RandomNURBSGenerator<PARAMETRIC_DIMENSIONALITY> &&other) = delete;
  RandomNURBSGenerator & operator=(const RandomNURBSGenerator<PARAMETRIC_DIMENSIONALITY> &rhs) = delete;
  RandomNURBSGenerator & operator=(RandomNURBSGenerator<PARAMETRIC_DIMENSIONALITY> &&rhs) = delete;
  ~RandomNURBSGenerator() override = default;

  RandomNURBSGenerator(std::array<ParametricCoordinate, 2> coord_limits, int max_degree, int dimension) {
    std::array<Degree, PARAMETRIC_DIMENSIONALITY> degrees =
        RandomSplineUtils<PARAMETRIC_DIMENSIONALITY>::GetRandomDegrees(max_degree);
    baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> knot_vectors =
        RandomSplineUtils<PARAMETRIC_DIMENSIONALITY>::GetRandomKnotVectors(coord_limits, degrees);
    std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points =
        RandomSplineUtils<PARAMETRIC_DIMENSIONALITY>::GetNumberOfPoints(degrees, knot_vectors);
    std::vector<src::spl::ControlPoint> control_points =
        RandomSplineUtils<PARAMETRIC_DIMENSIONALITY>::GetRandomControlPoints(dimension, number_of_points);
    std::vector<double> weights = RandomSplineUtils<PARAMETRIC_DIMENSIONALITY>::GetRandomWeights(number_of_points);
    this->physical_space_ = std::make_shared<src::spl::WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>>(
        control_points, weights, number_of_points);
    this->parameter_space_ =
        std::make_shared<src::spl::ParameterSpace<PARAMETRIC_DIMENSIONALITY>>(knot_vectors, degrees);
  }
};
}  // namespace splinelib::test::spl::random

#endif  // TEST_SPL_RANDOM_RANDOM_NURBS_GENERATOR_H_
