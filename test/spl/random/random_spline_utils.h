/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef TEST_SPL_RANDOM_RANDOM_SPLINE_UTILS_H_
#define TEST_SPL_RANDOM_RANDOM_SPLINE_UTILS_H_

#include <memory>
#include <vector>

#include "src/spl/b_spline.h"
#include "src/spl/nurbs.h"
#include "src/util/random.h"
#include "src/spl/parameter_space.h"
#include "src/spl/physical_space.h"

namespace splinelib::test {
using namespace splinelib::src;
template<int PARAMETRIC_DIMENSIONALITY>
class RandomSplineUtils {
 public:
  static std::shared_ptr<spl::BSpline<PARAMETRIC_DIMENSIONALITY>> GenerateRandomBSpline(
      std::array<ParametricCoordinate, 2> coord_limits, int max_degree, int dimension) {
    std::array<Degree, PARAMETRIC_DIMENSIONALITY> degrees = GetRandomDegrees(max_degree);
    baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> knot_vectors = GetRandomKnotVectors(coord_limits, degrees);
    std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points = GetNumberOfPoints(degrees, knot_vectors);
    std::vector<spl::ControlPoint> control_points = GetRandomControlPoints(dimension, number_of_points);
    auto physical_space = std::make_shared<spl::PhysicalSpace<PARAMETRIC_DIMENSIONALITY>>(control_points,
                                                                                          number_of_points);
    auto parameter_space = std::make_shared<spl::ParameterSpace<PARAMETRIC_DIMENSIONALITY>>(knot_vectors, degrees);
    return std::make_shared<spl::BSpline<PARAMETRIC_DIMENSIONALITY>>(physical_space, parameter_space);
  }

  static std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>> GenerateRandomNURBS(
      std::array<ParametricCoordinate, 2> coord_limits, int max_degree, int dimension) {
    std::array<Degree, PARAMETRIC_DIMENSIONALITY> degrees = GetRandomDegrees(max_degree);
    baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> knot_vectors = GetRandomKnotVectors(coord_limits, degrees);
    std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points = GetNumberOfPoints(degrees, knot_vectors);
    std::vector<src::spl::ControlPoint> control_points = GetRandomControlPoints(dimension, number_of_points);
    std::vector<double> weights = GetRandomWeights(number_of_points);
    auto physical_space = std::make_shared<src::spl::WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>>(
        control_points, weights, number_of_points);
    auto parameter_space = std::make_shared<src::spl::ParameterSpace<PARAMETRIC_DIMENSIONALITY>>(knot_vectors, degrees);
    return std::make_shared<spl::NURBS<PARAMETRIC_DIMENSIONALITY>>(physical_space, parameter_space);
  }

 private:
  static std::array<Degree, PARAMETRIC_DIMENSIONALITY> GetRandomDegrees(int max_degree) {
    std::array<Degree, PARAMETRIC_DIMENSIONALITY> degrees{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      degrees[i] = Degree{util::random::GetBinomialRandom<int>(1, max_degree, 1)};
    }
    return degrees;
  }

  static baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> GetRandomKnotVectors(
      std::array<ParametricCoordinate, 2> coord_limits, const std::array<Degree, PARAMETRIC_DIMENSIONALITY> &degree) {
    baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> knot_vectors;
    std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_knots = GetNumberOfKnots(degree);
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      std::vector<ParametricCoordinate> param_coord_vector;
      param_coord_vector.reserve(number_of_knots[i] - 2);
      for (int j = 0; j < degree[i].Get() + 1; ++j) {
        param_coord_vector.emplace_back(coord_limits[0]);
      }
      int number = number_of_knots[i] - 2 * degree[i].Get() - 2;
      for (int j = 1; j <= number; ++j) {
        double coord = (coord_limits[1].Get() - coord_limits[0].Get()) / (number + 1) * j + coord_limits[0].Get();
        param_coord_vector.emplace_back(coord);
      }
      for (int j = 0; j < degree[i].Get() + 1; ++j) {
        param_coord_vector.emplace_back(coord_limits[1]);
      }
      knot_vectors[i] = std::make_shared<baf::KnotVector>(param_coord_vector);
    }
    return knot_vectors;
  }

  static std::vector<src::spl::ControlPoint> GetRandomControlPoints(
      int dimension, const std::array<int, PARAMETRIC_DIMENSIONALITY> &number_of_points) {
    std::vector<src::spl::ControlPoint> control_points;
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(number_of_points);
    for (int i = 0; i < point_handler.GetNumberOfTotalMultiIndices(); ++i, ++point_handler) {
      std::vector<double> coordinates;
      coordinates.reserve(dimension);
      for (int j = 0; j < dimension; ++j) {
        coordinates.emplace_back(util::random::GetBinomialRandom<double>(-5, 5, 0.01));
      }
      control_points.emplace_back(coordinates);
    }
    return control_points;
  }

  static std::vector<double> GetRandomWeights(const std::array<int, PARAMETRIC_DIMENSIONALITY> &number_of_points) {
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(number_of_points);
    std::vector<double> weights;
    weights.reserve(point_handler.GetNumberOfTotalMultiIndices());
    for (int i = 0; i < point_handler.GetNumberOfTotalMultiIndices(); ++i, ++point_handler) {
      weights.emplace_back(util::random::GetBinomialRandom<double>(0.1, 2, 0.1));
    }
    return weights;
  }

  static std::array<int, PARAMETRIC_DIMENSIONALITY> GetNumberOfPoints(
      const std::array<Degree, PARAMETRIC_DIMENSIONALITY> &degrees,
      const baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> &knot_vectors) {
    std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      number_of_points[i] = knot_vectors[i]->GetNumberOfKnots() - degrees[i].Get() - 1;
    }
    return number_of_points;
  }

  static std::array<int, PARAMETRIC_DIMENSIONALITY> GetNumberOfKnots(
      const std::array<Degree, PARAMETRIC_DIMENSIONALITY> &degree) {
    std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_knots{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      number_of_knots[i] = util::random::GetBinomialRandom<int>(2 * degree[i].Get() + 2, 4 * degree[i].Get(), 1);
    }
    return number_of_knots;
  }
};
}  // namespace splinelib::test::spl::random

#endif  // TEST_SPL_RANDOM_RANDOM_SPLINE_UTILS_H_
