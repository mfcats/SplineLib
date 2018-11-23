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
#include "random.h"

namespace spl {
template<int DIM>
class RandomNURBSGenerator : public NURBSGenerator<DIM> {
 public:
  RandomNURBSGenerator() = default;
  virtual ~RandomNURBSGenerator() = default;

  RandomNURBSGenerator(std::array<ParamCoord, 2> param_coord_limits, int max_degree, int dimension) {
    std::array<Degree, DIM> degrees = GetRandomDegrees(max_degree);
    KnotVectors<DIM> knot_vectors = GetRandomKnotVectors(param_coord_limits, degrees);
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = knot_vectors[i]->GetNumberOfKnots() - degrees[i].get() - 1;
    }
    std::vector<baf::ControlPoint> control_points = GetRandomControlPoints(dimension, number_of_points);
    std::vector<double> weights = GetRandomWeights(number_of_points);
    this->physical_space_ = std::make_shared<WeightedPhysicalSpace<DIM>>(control_points, weights, number_of_points);
    this->parameter_space_ = std::make_shared<ParameterSpace<DIM>>(knot_vectors, degrees);
  }

 private:
  std::array<Degree, DIM> GetRandomDegrees(int max_degree) const {
    std::default_random_engine generator;
    std::binomial_distribution<int> distribution(max_degree, 0.5);
    std::array<Degree, DIM> degrees;
    for (int i = 0; i < DIM; ++i) {
      degrees[i] = Degree{util::Random::GetBinomialRandom<int>(0, max_degree, 1)};
    }
    return degrees;
  }

  KnotVectors<DIM> GetRandomKnotVectors(std::array<ParamCoord, 2> param_coord_limits,
                                        const std::array<Degree, DIM> &degrees) const {
    KnotVectors<DIM> knot_vectors;
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = util::Random::GetBinomialRandom<int>(2 * degrees[i].get() + 2, 4 * degrees[i].get(), 1);
    }
    std::array<std::vector<ParamCoord>, DIM> param_coord_vectors;
    for (int i = 0; i < DIM; ++i) {
      for (int j = 0; j < degrees[i].get() + 1; ++j) {
        param_coord_vectors[i].push_back(param_coord_limits[0]);
      }
      for (int j = 0; j < number_of_points[i] - 2 * degrees[i].get() - 2; ++j) {
        param_coord_vectors[i].push_back(ParamCoord((param_coord_limits[1].get() - param_coord_limits[1].get()) / 2));
      }
      for (int j = 0; j < degrees[i].get() + 1; ++j) {
        param_coord_vectors[i].push_back(param_coord_limits[1]);
      }
      baf::KnotVector knt_vec(param_coord_vectors[i]);
      knot_vectors[i] = std::make_shared<baf::KnotVector>(knt_vec);
    }
    return knot_vectors;
  }

  std::vector<baf::ControlPoint> GetRandomControlPoints(int dimension,
                                                        const std::array<int, DIM> &number_of_points) const {
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

  std::vector<double> GetRandomWeights(const std::array<int, DIM> &number_of_points) const {
    std::vector<double> weights;
    util::MultiIndexHandler<DIM> point_handler(number_of_points);
    for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
      weights.emplace_back(util::Random::GetBinomialRandom<double>(0.1, 5, 0.1));
    }
    return weights;
  }
};
}  // namespace spl

#endif  // SRC_SPL_RANDOM_NURBS_GENERATOR_H_
