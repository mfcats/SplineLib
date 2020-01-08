/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

// NOLINT(build/header_guard)

#include <vector>

template<int PARAMETRIC_DIMENSIONALITY>
PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::PhysicalSpace(std::vector<spl::ControlPoint> const &control_points,
                                                        std::array<int, PARAMETRIC_DIMENSIONALITY> const &number_of_points)
    : dimensionality_(control_points[0].GetDimensionality()), number_of_points_per_dimension_(number_of_points) {
  uint64_t total_number_of_points = 1;
  for (Dimension dimension{0}; dimension < Dimension{PARAMETRIC_DIMENSIONALITY}; ++dimension) {
    total_number_of_points *= number_of_points[dimension.Get()];
  }
  if (total_number_of_points != control_points.size()) {
    throw std::runtime_error(
        "The given number of control points in each dimension doesn't fit the length of the control point vector.");
  }
  for (auto &&control_point : control_points) {
    if (control_point.GetDimensionality() != dimensionality_) {
      throw std::runtime_error("The dimension has to be the same for all control points.");
    }
    for (Dimension dimension{0}; dimension < Dimension{dimensionality_}; ++dimension) {
      control_points_.emplace_back(control_point[dimension]);
    }
  }
}

template<int PARAMETRIC_DIMENSIONALITY>
bool PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AreEqual(PhysicalSpace<PARAMETRIC_DIMENSIONALITY> const &rhs,
                                                        Tolerance const &tolerance) const {
  return std::equal(
      control_points_.begin(), control_points_.end(), rhs.control_points_.begin(), rhs.control_points_.end(),
      [&](double cp_a, double cp_b) {return util::numeric_settings::AreEqual<double>(cp_a, cp_b, tolerance.Get());}) &&
         std::equal(
             number_of_points_per_dimension_.begin(), number_of_points_per_dimension_.end(),
             rhs.number_of_points_per_dimension_.begin(), rhs.number_of_points_per_dimension_.end(),
             [&](int number_a, int number_b) {return util::numeric_settings::AreEqual<double>(number_a, number_b);});
}

template<int PARAMETRIC_DIMENSIONALITY>
spl::ControlPoint PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetControlPoint(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) const {
  auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(multi_index);
  int first = dimensionality_ * point_handler.GetCurrent1DIndex();
  std::vector<double> coordinates;
  coordinates.reserve(dimensionality_);
  for (int coordinate = 0; coordinate < dimensionality_; coordinate++) {
    coordinates.emplace_back(control_points_[first + coordinate]);
  }
  return spl::ControlPoint(coordinates);
}

template<int PARAMETRIC_DIMENSIONALITY>
spl::ControlPoint PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetControlPoint(int index_1d) const {
  auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(index_1d);
  int first = dimensionality_ * point_handler.GetCurrent1DIndex();
  std::vector<double> coordinates;
  coordinates.reserve(dimensionality_);
  for (int coordinate = 0; coordinate < dimensionality_; coordinate++) {
    coordinates.emplace_back(control_points_[first + coordinate]);
  }
  return spl::ControlPoint(coordinates);
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetControlPoint(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index,
    ControlPoint const &control_point, Dimension const &dimension, int (*before)(int)) {
  const std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points_before(number_of_points_per_dimension_);
  if (before != nullptr)
    GetValue(number_of_points_per_dimension_, dimension) = before(GetValue(number_of_points_per_dimension_, dimension));
  auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(multi_index);
  int first = dimensionality_ * point_handler.GetCurrent1DIndex();
  for (int coordinate = 0; coordinate < dimensionality_; coordinate++) {
    control_points_[first + coordinate] = control_point[Dimension{coordinate}];
  }
  number_of_points_per_dimension_ = number_of_points_before;
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetControlPoint(int index_1d, ControlPoint const &control_point,
                                                               Dimension const &dimension, int (*before)(int)) {
  const std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points_before(number_of_points_per_dimension_);
  if (before != nullptr)
    GetValue(number_of_points_per_dimension_, dimension) = before(GetValue(number_of_points_per_dimension_, dimension));
  auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(index_1d);
  int first = dimensionality_ * point_handler.GetCurrent1DIndex();
  for (int coordinate = 0; coordinate < dimensionality_; coordinate++) {
    control_points_[first + coordinate] = control_point[Dimension{coordinate}];
  }
  number_of_points_per_dimension_ = number_of_points_before;
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AddControlPoints(int number) {
  for (int i = 0; i < number; ++i) {
    for (int j = 0; j < dimensionality_; ++j) {
      control_points_.emplace_back(0.0);
    }
  }
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::RemoveControlPoints(int number) {
  control_points_.erase(control_points_.end() - number * dimensionality_, control_points_.end());
}

template<int PARAMETRIC_DIMENSIONALITY>
int PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetTotalNumberOfControlPoints() const {
  return static_cast<int>(control_points_.size()) / dimensionality_;
}

template<int PARAMETRIC_DIMENSIONALITY>
std::array<int, PARAMETRIC_DIMENSIONALITY> PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetNumberOfPointsPerDirection() const {
  return number_of_points_per_dimension_;
}

template<int PARAMETRIC_DIMENSIONALITY>
std::array<int, PARAMETRIC_DIMENSIONALITY>
PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetMaximumPointIndexPerDirection() const {
  std::array<int, PARAMETRIC_DIMENSIONALITY> maximum_index = number_of_points_per_dimension_;
  for (auto &index : maximum_index) {
    --index;
  }
  return maximum_index;
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::IncrementNumberOfPoints(int dimension) {
  ++number_of_points_per_dimension_[dimension];
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::DecrementNumberOfPoints(int dimension) {
  --number_of_points_per_dimension_[dimension];
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetNumberOfPoints(int dimension, int value) {
  number_of_points_per_dimension_[dimension] = value;
}

template<int PARAMETRIC_DIMENSIONALITY>
int PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetDimensionality() const {
  return dimensionality_;
}

template<int PARAMETRIC_DIMENSIONALITY>
Weight PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetWeight(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &/*indices*/) const {
  return Weight{1.0};
}

template<int PARAMETRIC_DIMENSIONALITY>
Weight PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetWeight(int /*index_1d*/) const {
  return Weight{1.0};
}

template<int PARAMETRIC_DIMENSIONALITY>
double PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetExpansion() const {
  double expansion = 0;
  for (auto &cp : control_points_) {
    if (cp > expansion) expansion = cp;
  }
  return expansion;
}

template<int PARAMETRIC_DIMENSIONALITY>
double PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetMaximumDistanceFromOrigin() const {
  double maximum = 0;
  std::vector<double> point(dimensionality_, 0);
  std::vector<double> origin(dimensionality_, 0);
  for (int i = 0; i < GetTotalNumberOfControlPoints(); ++i) {
    for (int j = 0; j < dimensionality_; ++j) {
      point[j] = control_points_[i * dimensionality_ + j];
    }
    double distance = util::vector_utils::ComputeDistance(origin, point);
    if (distance > maximum) maximum = distance;
  }
  return maximum;
}

template<int PARAMETRIC_DIMENSIONALITY>
std::vector<spl::ControlPoint> PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SplitControlPoints(const Dimension &dimension,
                                                                                            int offset,
                                                                                            int length) {
  std::vector<spl::ControlPoint> points;
  std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = GetNumberOfPointsPerDirection();
  GetValue(point_handler_length, dimension) = length;
  util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
  for (int i = 0; i < point_handler.GetNumberOfTotalMultiIndices(); ++i, ++point_handler) {
    auto indices = point_handler.GetCurrentIndex();
    GetValue(indices, dimension) += offset;
    points.emplace_back(GetControlPoint(indices));
  }
  return points;
}

template<int PARAMETRIC_DIMENSIONALITY>
int PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetNumberOfPointsForDimension(Dimension const &dimension) const {
  return GetValue(number_of_points_per_dimension_, dimension);
}

template<int PARAMETRIC_DIMENSIONALITY>
std::vector<double> PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetControlPoints() const {
  return control_points_;
}
