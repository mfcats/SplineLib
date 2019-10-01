/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_RANDOM_SPLINE_UTILS_H_
#define SRC_SPL_RANDOM_SPLINE_UTILS_H_

#include <memory>
#include <vector>

#include "random.h"
#include "parameter_space.h"
#include "physical_space.h"

namespace spl {
template<int DIM>
class RandomSplineUtils {
 public:
  RandomSplineUtils() = default;
  virtual ~RandomSplineUtils() = default;

  static std::array<Degree, DIM> GetRandomDegrees(int max_degree) {
    std::array<Degree, DIM> degrees;
    for (int i = 0; i < DIM; ++i) {
      degrees[i] = Degree{util::Random::GetBinomialRandom<int>(1, max_degree, 1)};
    }
    return degrees;
  }

  static KnotVectors<DIM> GetRandomKnotVectors(std::array<ParamCoord, 2> coord_limits,
                                               const std::array<Degree, DIM> &degree) {
    KnotVectors<DIM> knot_vectors;
    std::array<int, DIM> number_of_knots = GetNumberOfKnots(degree);
    for (int i = 0; i < DIM; ++i) {
      std::vector<ParamCoord> param_coord_vector;
      for (int j = 0; j < degree[i].get() + 1; ++j) {
        param_coord_vector.push_back(coord_limits[0]);
      }
      int number = number_of_knots[i] - 2 * degree[i].get() - 2;
      for (int j = 1; j <= number; ++j) {
        double coord = (coord_limits[1].get() - coord_limits[0].get()) / (number + 1) * j + coord_limits[0].get();
        param_coord_vector.emplace_back(coord);
      }
      for (int j = 0; j < degree[i].get() + 1; ++j) {
        param_coord_vector.push_back(coord_limits[1]);
      }
      knot_vectors[i] = std::make_shared<baf::KnotVector>(param_coord_vector);
    }
    return knot_vectors;
  }

  static std::vector<baf::ControlPoint> GetRandomControlPoints(int dimension,
                                                               const std::array<int, DIM> &number_of_points) {
    std::vector<baf::ControlPoint> control_points;
    util::MultiIndexHandler<DIM> point_handler(number_of_points);
    for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
      std::vector<double> coordinates;
      for (int j = 0; j < dimension; ++j) {
        coordinates.push_back(util::Random::GetBinomialRandom<double>(-5, 5, 0.01));
      }
      control_points.emplace_back(coordinates);
    }
    return control_points;
  }

  static std::vector<double> GetRandomWeights(const std::array<int, DIM> &number_of_points) {
    std::vector<double> weights;
    util::MultiIndexHandler<DIM> point_handler(number_of_points);
    for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
      weights.emplace_back(util::Random::GetBinomialRandom<double>(0.1, 2, 0.1));
    }
    return weights;
  }

  static std::array<int, DIM> GetNumberOfPoints(const std::array<Degree, DIM> &degrees,
                                                const KnotVectors<DIM> &knot_vectors) {
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = knot_vectors[i]->GetNumberOfKnots() - degrees[i].get() - 1;
    }
    return number_of_points;
  }

 private:
  static std::array<int, DIM> GetNumberOfKnots(const std::array<Degree, DIM> &degree) {
    std::array<int, DIM> number_of_knots;
    for (int i = 0; i < DIM; ++i) {
      number_of_knots[i] = util::Random::GetBinomialRandom<int>(2 * degree[i].get() + 2, 4 * degree[i].get(), 1);
    }
    return number_of_knots;
  }
};
}  // namespace spl

#endif  // SRC_SPL_RANDOM_SPLINE_UTILS_H_
