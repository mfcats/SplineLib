/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include <vector>


template<int PARAMETRIC_DIMENSIONALITY>
WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::WeightedPhysicalSpace(
    std::vector<spl::ControlPoint> const &control_points, std::vector<Weight> const &weights,
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &number_of_points)
    : PhysicalSpace<PARAMETRIC_DIMENSIONALITY>(control_points, number_of_points), weights_(weights) {
  if (control_points.size() != weights_.size()) {
    throw std::runtime_error("splinelib::src::spl::WeightedPhysicalSpace::WeightedPhysicalSpace: The number of control"
                             "points and weights have to be the same.");
  }
}

template<int PARAMETRIC_DIMENSIONALITY>
Weight WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetWeight(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) const {
  util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(this->number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(multi_index);
  return GetValue(weights_, point_handler.GetCurrent1DIndex());
}

template<int PARAMETRIC_DIMENSIONALITY>
Weight WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetWeight(int index_1d) const {
  return GetValue(weights_, index_1d);
}

template<int PARAMETRIC_DIMENSIONALITY>
std::vector<Weight> WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetWeights() const {
  return weights_;
}

template<int PARAMETRIC_DIMENSIONALITY>
void WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetWeight(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index, Weight const &weight, Dimension const &dimension,
    int (*before)(int)) {
  const std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points_before(this->number_of_points_per_dimension_);
  if (before != nullptr) {
    GetValue(this->number_of_points_per_dimension_, dimension) =
        before(GetValue(this->number_of_points_per_dimension_, dimension));
  }
  util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(this->number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(multi_index);
  GetValue(weights_, point_handler.GetCurrent1DIndex()) = weight;
  this->number_of_points_per_dimension_ = number_of_points_before;
}

template<int PARAMETRIC_DIMENSIONALITY>
void WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetWeight(int index_1d, Weight const &weight,
                                                                 Dimension const &dimension, int (*before)(int)) {
  const std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points_before(this->number_of_points_per_dimension_);
  if (before != nullptr) {
    GetValue(this->number_of_points_per_dimension_, dimension) =
        before(GetValue(this->number_of_points_per_dimension_, dimension));
  }
  GetValue(weights_, index_1d) = weight;
  this->number_of_points_per_dimension_ = number_of_points_before;
}

template<int PARAMETRIC_DIMENSIONALITY>
ControlPoint WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetHomogenousControlPoint(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) const {
  std::vector<double> coordinates(this->dimensionality_ + 1, 0.0);
  util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(this->number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(multi_index);
  int index_first_coordinate = this->dimensionality_ * point_handler.GetCurrent1DIndex();
  double const weight = GetValue(weights_, point_handler.GetCurrent1DIndex()).Get();
  for (Dimension current_dimension{0}; current_dimension < Dimension{this->dimensionality_}; ++current_dimension) {
    int index_1d_current_coordinate = index_first_coordinate + current_dimension.Get();
    double homogeneous_coordinate = GetValue(this->control_points_, index_1d_current_coordinate) * weight;
    GetValue(coordinates, current_dimension) = homogeneous_coordinate;
  }
  GetValue(coordinates, this->dimensionality_) = weight;
  return ControlPoint(coordinates);
}

template<int PARAMETRIC_DIMENSIONALITY>
ControlPoint WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetHomogenousControlPoint(int index_1d) const {
  std::vector<double> coordinates(this->dimensionality_ + 1, 0.0);
  int index_first_coordinate = this->dimensionality_ * index_1d;
  double const weight = GetValue(weights_, index_1d).Get();
  for (Dimension current_dimension{0}; current_dimension < Dimension{this->dimensionality_}; ++current_dimension) {
    int index_1d_current_coordinate = index_first_coordinate + current_dimension.Get();
    double homogeneous_coordinate = GetValue(this->control_points_, index_1d_current_coordinate) * weight;
    GetValue(coordinates, current_dimension) = homogeneous_coordinate;
  }
  GetValue(coordinates, this->dimensionality_) = weight;
  return ControlPoint(coordinates);
}

template<int PARAMETRIC_DIMENSIONALITY>
bool WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AreEqual(WeightedPhysicalSpace const &rhs,
                                                                Tolerance const &tolerance) const {
  return std::equal(weights_.begin(), weights_.end(), rhs.weights_.begin(), rhs.weights_.end(),
                    [&](Weight const &weight_a, Weight const &weight_b) {
                      return util::numeric_settings::AreEqual<double>(weight_a.Get(), weight_b.Get(), tolerance.Get());
                    }) &&
         PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AreEqual(rhs, Tolerance{tolerance});
}

template<int PARAMETRIC_DIMENSIONALITY>
void WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetWeightedControlPoint(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index,
    ControlPoint const &control_point, Weight const &weight) {
  PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetControlPoint(multi_index, control_point);
  SetWeight(multi_index, weight);
}

template<int PARAMETRIC_DIMENSIONALITY>
double WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetMinimumWeight() const {
  double minimum_weight_value = GetFront(weights_).Get();
  for (const auto &weight : weights_) {
    if (weight.Get() < minimum_weight_value) minimum_weight_value = weight.Get();
  }
  return minimum_weight_value;
}

template<int PARAMETRIC_DIMENSIONALITY>
std::vector<Weight> WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SplitWeights(Dimension dimension,
                                                                                   int offset,
                                                                                   int length) {
  std::array<int, PARAMETRIC_DIMENSIONALITY> weight_handler_length = this->GetNumberOfPointsPerDirection();
  GetValue(weight_handler_length, dimension) = length;
  util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> weight_handler(weight_handler_length);
  std::vector<Weight> weights(weight_handler.GetNumberOfTotalMultiIndices(), Weight{0.0});
  for (int weight_index = 0; weight_index < weight_handler.GetNumberOfTotalMultiIndices();
       ++weight_index, ++weight_handler) {
    auto multi_index = weight_handler.GetCurrentIndex();
    GetValue(multi_index, dimension) += offset;
    GetValue(weights, weight_index) = GetWeight(multi_index);
  }
  return weights;
}

template<int PARAMETRIC_DIMENSIONALITY>
void WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AddControlPoints(int number) {
  PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AddControlPoints(number);
  for (int i = 0; i < number; ++i) weights_.emplace_back(Weight{0.0});
}

template<int PARAMETRIC_DIMENSIONALITY>
void WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::RemoveControlPoints(int number) {
  PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::RemoveControlPoints(number);
  auto last = weights_.end();
  weights_.erase(last - number, last);
}
