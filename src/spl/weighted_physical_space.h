/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_SPL_WEIGHTED_PHYSICAL_SPACE_H_
#define SRC_SPL_WEIGHTED_PHYSICAL_SPACE_H_

#include <utility>
#include <vector>

#include "src/spl/physical_space.h"

namespace splinelib::src::spl {
template<int PARAMETRIC_DIMENSIONALITY>
class WeightedPhysicalSpace : public PhysicalSpace<PARAMETRIC_DIMENSIONALITY> {
 public:
  WeightedPhysicalSpace() = default;
  WeightedPhysicalSpace(const std::vector<spl::ControlPoint> &control_points, std::vector<double> weights,
      std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points)
      : PhysicalSpace<PARAMETRIC_DIMENSIONALITY>(control_points, number_of_points), weights_(std::move(weights)) {
    if (control_points.size() != weights_.size()) {
      throw std::runtime_error("The number of control points and weights has to be the same.");
    }
  }

  WeightedPhysicalSpace(const WeightedPhysicalSpace &physical_space)
      : PhysicalSpace<PARAMETRIC_DIMENSIONALITY>(physical_space) {
    this->weights_ = physical_space.weights_;
  }

  WeightedPhysicalSpace(WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY> &&other) noexcept = default;
  WeightedPhysicalSpace & operator=(const WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY> &rhs) = default;
  WeightedPhysicalSpace & operator=(WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY> &&rhs) noexcept = default;
  ~WeightedPhysicalSpace() override = default;

  bool AreEqual(const WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY> &rhs,
                double tolerance = util::numeric_settings::GetEpsilon<double>()) const {
    return std::equal(weights_.begin(), weights_.end(),
                      rhs.weights_.begin(), rhs.weights_.end(),
                      [&](double weight_a, double weight_b) {
                        return util::numeric_settings::AreEqual<double>(weight_a, weight_b, tolerance);
                      }) && PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AreEqual(rhs, tolerance);
  }

  virtual spl::ControlPoint GetHomogenousControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> indices) const {
    auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(this->number_of_points_);
    point_handler.SetCurrentIndex(indices);
    int first = this->dimensionality_ * point_handler.GetCurrent1DIndex();
    std::vector<double> coordinates;
    coordinates.reserve(this->dimensionality_);
    for (int coordinate = 0; coordinate < this->dimensionality_; coordinate++) {
      coordinates.emplace_back(this->control_points_[first + coordinate] * weights_[first / this->dimensionality_]);
    }
    coordinates.push_back(this->weights_[point_handler.GetCurrent1DIndex()]);
    return spl::ControlPoint(coordinates);
  }

  virtual spl::ControlPoint GetHomogenousControlPoint(int index_1d) const {
    auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(this->number_of_points_);
    point_handler.SetCurrentIndex(index_1d);
    int first = this->dimensionality_ * point_handler.GetCurrent1DIndex();
    std::vector<double> coordinates;
    coordinates.reserve(this->dimensionality_);
    for (int coordinate = 0; coordinate < this->dimensionality_; coordinate++) {
      coordinates.emplace_back(this->control_points_[first + coordinate] * weights_[first / this->dimensionality_]);
    }
    coordinates.push_back(this->weights_[point_handler.GetCurrent1DIndex()]);
    return spl::ControlPoint(coordinates);
  }

  double GetWeight(std::array<int, PARAMETRIC_DIMENSIONALITY> indices) const override {
    auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(this->number_of_points_);
    point_handler.SetCurrentIndex(indices);
    int first = point_handler.GetCurrent1DIndex();
    return weights_[first];
  }

  double GetWeight(int index_1d) const override {
    auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(this->number_of_points_);
    point_handler.SetCurrentIndex(index_1d);
    int first = point_handler.GetCurrent1DIndex();
    return weights_[first];
  }

  double GetMinimumWeight() const {
    double minimum = weights_[0];
    for (const auto &weight : weights_) {
      if (weight < minimum) minimum = weight;
    }
    return minimum;
  }

  void SetWeightedControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> indices,
      const spl::ControlPoint &control_point, double weight) {
    PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetControlPoint(indices, control_point);
    SetWeight(indices, weight);
  }

  void SetWeight(std::array<int, PARAMETRIC_DIMENSIONALITY> indices, double weight, int dimension = 0,
      int (*before)(int) = nullptr) {
    const std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points_before(this->number_of_points_);
    if (before != nullptr) this->number_of_points_[dimension] = before(this->number_of_points_[dimension]);
    auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(this->number_of_points_);
    point_handler.SetCurrentIndex(indices);
    int first = point_handler.GetCurrent1DIndex();
    weights_[first] = weight;
    this->number_of_points_ = number_of_points_before;
  }

  void SetWeight(int index_1d, double weight, int dimension = 0,
                 int (*before)(int) = nullptr) {
    const std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points_before(this->number_of_points_);
    if (before != nullptr) this->number_of_points_[dimension] = before(this->number_of_points_[dimension]);
    auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(this->number_of_points_);
    point_handler.SetCurrentIndex(index_1d);
    int first = point_handler.GetCurrent1DIndex();
    weights_[first] = weight;
    this->number_of_points_ = number_of_points_before;
  }

  void AddControlPoints(int number) override {
    PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AddControlPoints(number);
    for (int i = 0; i < number; ++i) {
      weights_.emplace_back(0.0);
    }
  }

  void RemoveControlPoints(int number) override {
    PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::RemoveControlPoints(number);
    weights_.erase(weights_.end() - number, weights_.end());
  }

  std::vector<double> GetWeights() const {
    return weights_;
  }

  std::vector<double> GetDividedWeights(int first, int length, int dimension) {
    std::vector<double> weights;
    std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = this->GetPointsPerDirection();
    point_handler_length[dimension] = length;
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
    for (int i = 0; i < point_handler.GetNumberOfTotalMultiIndices(); ++i, ++point_handler) {
      auto indices = point_handler.GetCurrentIndex();
      indices[dimension] += first;
      weights.emplace_back(GetWeight(indices));
    }
    return weights;
  }

 private:
  std::vector<double> weights_;
};
}  // namespace splinelib::src::spl

#endif  // SRC_SPL_WEIGHTED_PHYSICAL_SPACE_H_
