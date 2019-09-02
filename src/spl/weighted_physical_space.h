/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_WEIGHTED_PHYSICAL_SPACE_H_
#define SRC_SPL_WEIGHTED_PHYSICAL_SPACE_H_

#include <vector>

#include "physical_space.h"

namespace spl {
template<int DIM>
class WeightedPhysicalSpace : public PhysicalSpace<DIM> {
 public:
  WeightedPhysicalSpace() = default;
  WeightedPhysicalSpace(const std::vector<baf::ControlPoint> &control_points,
                        const std::vector<double> &weights,
                        std::array<int, DIM> number_of_points) : PhysicalSpace<DIM>(control_points,
                                                                                    number_of_points),
                                                                 weights_(weights) {
    if (control_points.size() != weights_.size()) {
      throw std::runtime_error("The number of control points and weights has to be the same.");
    }
  }

  WeightedPhysicalSpace(const WeightedPhysicalSpace &physical_space) : PhysicalSpace<DIM>(physical_space) {
    this->weights_ = physical_space.weights_;
  }

  bool AreEqual(const WeightedPhysicalSpace<DIM> &rhs,
                double tolerance = util::NumericSettings<double>::kEpsilon()) const {
    return std::equal(weights_.begin(), weights_.end(),
                      rhs.weights_.begin(), rhs.weights_.end(),
                      [&](double weight_a, double weight_b) {
                        return util::NumericSettings<double>::AreEqual(weight_a, weight_b, tolerance);
                      }) && PhysicalSpace<DIM>::AreEqual(rhs, tolerance);
  }

  virtual baf::ControlPoint GetHomogenousControlPoint(std::array<int, DIM> indices) const {
    std::vector<double> coordinates;
    util::MultiIndexHandler<DIM> point_handler = util::MultiIndexHandler<DIM>(this->number_of_points_);
    point_handler.SetIndices(indices);
    int first = this->dimension_ * point_handler.Get1DIndex();
    for (int coordinate = 0; coordinate < this->dimension_; coordinate++) {
      coordinates.push_back(this->control_points_[first + coordinate] * weights_[first / this->dimension_]);
    }
    coordinates.push_back(this->weights_[point_handler.Get1DIndex()]);
    return baf::ControlPoint(coordinates);
  }

  double GetWeight(std::array<int, DIM> indices) const override {
    util::MultiIndexHandler<DIM> point_handler = util::MultiIndexHandler<DIM>(this->number_of_points_);
    point_handler.SetIndices(indices);
    int first = point_handler.Get1DIndex();
    return weights_[first];
  }

  double GetMinimumWeight() const {
    double minimum = weights_[0];
    for (const auto &weight : weights_) {
      if (weight < minimum) minimum = weight;
    }
    return minimum;
  }

  void SetWeightedControlPoint(std::array<int, DIM> indices, const baf::ControlPoint &control_point, double weight) {
    PhysicalSpace<DIM>::SetControlPoint(indices, control_point);
    SetWeight(indices, weight);
  }

  void SetWeight(std::array<int, DIM> indices, double weight, int dimension = 0, int (*before)(int) = nullptr) {
    const std::array<int, DIM> number_of_points_before(this->number_of_points_);
    if (before) this->number_of_points_[dimension] = before(this->number_of_points_[dimension]);
    util::MultiIndexHandler<DIM> point_handler = util::MultiIndexHandler<DIM>(this->number_of_points_);
    point_handler.SetIndices(indices);
    int first = point_handler.Get1DIndex();
    weights_[first] = weight;
    this->number_of_points_ = number_of_points_before;
  }

  void AddControlPoints(int number) override {
    PhysicalSpace<DIM>::AddControlPoints(number);
    for (int i = 0; i < number; ++i) {
      weights_.emplace_back(0.0);
    }
  }

  void RemoveControlPoints(int number) override {
    PhysicalSpace<DIM>::RemoveControlPoints(number);
    weights_.erase(weights_.end() - number, weights_.end());
  }

  std::vector<double> GetWeights() const {
    return weights_;
  }

  std::vector<double> GetDividedWeights(int first, int length, int dimension) {
    std::vector<double> weights;
    std::array<int, DIM> point_handler_length = this->GetPointsPerDirection();
    point_handler_length[dimension] = length;
    util::MultiIndexHandler<DIM> point_handler(point_handler_length);
    for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
      auto indices = point_handler.GetIndices();
      indices[dimension] += first;
      weights.emplace_back(GetWeight(indices));
    }
    return weights;
  }

 private:
  std::vector<double> weights_;
};
}  // namespace spl

#endif  // SRC_SPL_WEIGHTED_PHYSICAL_SPACE_H_
