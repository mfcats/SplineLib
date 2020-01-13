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
PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::PhysicalSpace(
    std::vector<ControlPoint> const &control_points,
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &number_of_points_per_direction)
    : dimensionality_(GetFront(control_points).GetDimensionality()),
      number_of_points_per_dimension_(number_of_points_per_direction) {
  total_number_of_points_ = 1;
  for (Dimension dimension{0}; dimension < Dimension{PARAMETRIC_DIMENSIONALITY}; ++dimension) {
    total_number_of_points_ *= GetValue(number_of_points_per_dimension_, dimension);
  }
  if (total_number_of_points_ != control_points.size()) {
    throw std::runtime_error(
        "splinelib::src::spl::PhysicalSpace::PhysicalSpace: Multiplying the specified number of points per dimension"
        "does not yield the total number of points provided.");
  }
  for (auto &&control_point : control_points) {
    if (control_point.GetDimensionality() != dimensionality_) {
      throw std::runtime_error("splinelib::src::spl::PhysicalSpace::PhysicalSpace: Control point of insufficient"
                               "dimensionality provided.");
    }
    for (Dimension dimension{0}; dimension < Dimension{dimensionality_}; ++dimension) {
      control_points_.emplace_back(control_point[dimension]);
    }
  }
}

template<int PARAMETRIC_DIMENSIONALITY>
int PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetDimensionality() const {
  return dimensionality_;
}

template<int PARAMETRIC_DIMENSIONALITY>
int PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetNumberOfPointsForDimension(Dimension const &dimension) const {
  return GetValue(number_of_points_per_dimension_, dimension);
}

template<int PARAMETRIC_DIMENSIONALITY>
std::array<int, PARAMETRIC_DIMENSIONALITY>
PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetNumberOfPointsPerDirection() const {
  return number_of_points_per_dimension_;
}

template<int PARAMETRIC_DIMENSIONALITY>
std::array<int, PARAMETRIC_DIMENSIONALITY>
PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetMaximumPointIndexPerDirection() const {
  std::array<int, PARAMETRIC_DIMENSIONALITY> maximum_index = number_of_points_per_dimension_;
  for (auto &index : maximum_index) --index;
  return maximum_index;
}

template<int PARAMETRIC_DIMENSIONALITY>
int PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetTotalNumberOfControlPoints() const {
  return static_cast<int>(control_points_.size()) / dimensionality_;
}

template<int PARAMETRIC_DIMENSIONALITY>
spl::ControlPoint PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetControlPoint(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) const {
  std::vector<double> coordinates(dimensionality_, 0.0);
  util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(multi_index);
  for (Dimension current_dimension{0}; current_dimension < Dimension{dimensionality_}; ++current_dimension) {
    int const index_1d_current_coordinate = (dimensionality_ * point_handler.GetCurrent1DIndex()) +
                                            current_dimension.Get();
    GetValue(coordinates, current_dimension) = GetValue(control_points_, index_1d_current_coordinate);
  }
  return ControlPoint(coordinates);
}

template<int PARAMETRIC_DIMENSIONALITY>
spl::ControlPoint PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetControlPoint(int index_1d) const {
  std::vector<double> coordinates(dimensionality_, 0.0);
  for (Dimension current_dimension{0}; current_dimension < Dimension{dimensionality_}; ++current_dimension) {
    int const index_1d_current_coordinate = (dimensionality_ * index_1d) + current_dimension.Get();
    GetValue(coordinates, current_dimension) = GetValue(control_points_, index_1d_current_coordinate);
  }
  return ControlPoint(coordinates);
}

template<int PARAMETRIC_DIMENSIONALITY>
std::vector<double> PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetControlPoints() const {
  return control_points_;
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
std::vector<Weight> PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetWeights() const {
  return std::vector<Weight>(total_number_of_points_, Weight{1.0});
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetControlPoint(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index, ControlPoint const &control_point,
    Dimension const &dimension, int (*before)(int)) {
  const std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points_before(number_of_points_per_dimension_);
  if (before != nullptr) {
    GetValue(number_of_points_per_dimension_, dimension) = before(GetValue(number_of_points_per_dimension_, dimension));
  }
  util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(multi_index);
  for (Dimension current_dimension{0}; current_dimension < Dimension{dimensionality_}; ++current_dimension) {
    int const index_1d_current_coordinate = (dimensionality_ * point_handler.GetCurrent1DIndex()) +
                                            current_dimension.Get();
    GetValue(control_points_, index_1d_current_coordinate) = control_point[current_dimension];
  }
  number_of_points_per_dimension_ = number_of_points_before;
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetControlPoint(int index_1d, ControlPoint const &control_point,
                                                               Dimension const &dimension, int (*before)(int)) {
  const std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points_before(number_of_points_per_dimension_);
  if (before != nullptr) {
    GetValue(number_of_points_per_dimension_, dimension) = before(GetValue(number_of_points_per_dimension_, dimension));
  }
  for (Dimension current_dimension{0}; current_dimension < Dimension{dimensionality_}; ++current_dimension) {
    int const index_1d_current_coordinate = (dimensionality_ * index_1d) + current_dimension.Get();
    GetValue(control_points_, index_1d_current_coordinate) = control_point[current_dimension];
  }
  number_of_points_per_dimension_ = number_of_points_before;
}

template<int PARAMETRIC_DIMENSIONALITY>
double PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetExpansion() const {
  double expansion = 0;
  for (auto &current_control_point : control_points_) {
    if (current_control_point > expansion) expansion = current_control_point;
  }
  return expansion;
}

template<int PARAMETRIC_DIMENSIONALITY>
double PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetMaximumDistanceFromOrigin() const {
  double maximum_distance_from_origin = 0;
  std::vector<double> current_point(dimensionality_, 0);
  std::vector<double> origin(dimensionality_, 0);
  for (int current_point_index = 0; current_point_index < GetTotalNumberOfControlPoints(); ++current_point_index) {
    for (Dimension current_dimension{0}; current_dimension < Dimension{dimensionality_}; ++current_dimension) {
      int const index_1d_current_coordinate = (dimensionality_ * current_point_index) + current_dimension.Get();
      GetValue(current_point, current_dimension) = GetValue(control_points_, index_1d_current_coordinate);
    }
    double current_distance = util::vector_utils::ComputeDistance(origin, current_point);
    if (current_distance > maximum_distance_from_origin) maximum_distance_from_origin = current_distance;
  }
  return maximum_distance_from_origin;
}

template<int PARAMETRIC_DIMENSIONALITY>
std::vector<spl::ControlPoint> PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SplitControlPoints(
    Dimension const &dimension, int offset, int length) {
  std::vector<spl::ControlPoint> points;
  std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = GetNumberOfPointsPerDirection();
  GetValue(point_handler_length, dimension) = length;
  util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
  for (int i = 0; i < point_handler.GetNumberOfTotalMultiIndices(); ++i, ++point_handler) {
    auto multi_index = point_handler.GetCurrentIndex();
    GetValue(multi_index, dimension) += offset;
    points.emplace_back(GetControlPoint(multi_index));
  }
  return points;
}

template<int PARAMETRIC_DIMENSIONALITY>
bool PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AreEqual(PhysicalSpace<PARAMETRIC_DIMENSIONALITY> const &rhs,
                                                        Tolerance const &tolerance) const {
  return std::equal(control_points_.begin(), control_points_.end(),
                    rhs.control_points_.begin(), rhs.control_points_.end(),
                    [&](double cp_a, double cp_b) {
                        return util::numeric_settings::AreEqual<double>(cp_a, cp_b, tolerance.Get());
                    }) &&
         std::equal(number_of_points_per_dimension_.begin(), number_of_points_per_dimension_.end(),
                    rhs.number_of_points_per_dimension_.begin(), rhs.number_of_points_per_dimension_.end(),
                    [&](int number_a, int number_b) {
                        return util::numeric_settings::AreEqual<double>(number_a, number_b);
                    });
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::IncrementNumberOfPoints(Dimension dimension) {
  ++GetValue(number_of_points_per_dimension_, dimension);
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::DecrementNumberOfPoints(Dimension dimension) {
  --GetValue(number_of_points_per_dimension_, dimension);
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AddControlPoints(int number) {
  for (int i = 0; i < number; ++i) {
    for (Dimension current_dimension{0}; current_dimension < Dimension{dimensionality_}; ++current_dimension) {
      control_points_.emplace_back(0.0);
    }
  }
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::RemoveControlPoints(int number) {
  auto last = control_points_.end();
  control_points_.erase(last - number * dimensionality_, last);
}

template<int PARAMETRIC_DIMENSIONALITY>
void PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetNumberOfPoints(Dimension dimension, int value) {
  GetValue(number_of_points_per_dimension_, dimension) = value;
}
