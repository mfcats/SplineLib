/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_PHYSICAL_SPACE_H_
#define SRC_SPL_PHYSICAL_SPACE_H_

#include <stdexcept>
#include <vector>

#include "control_point.h"
#include "multi_index_handler.h"

namespace spl {
template<int DIM>
class PhysicalSpace {
 public:
  PhysicalSpace() = default;
  virtual ~PhysicalSpace() = default;
  explicit PhysicalSpace(const std::vector<baf::ControlPoint> &control_points, std::array<int, DIM> number_of_points)
      : dimension_(control_points[0].GetDimension()), number_of_points_(number_of_points) {
    uint64_t total_number_of_points = 1;
    for (int dim = 0; dim < DIM; dim++) {
      total_number_of_points *= number_of_points[dim];
    }
    if (total_number_of_points != control_points.size()) {
      throw std::runtime_error(
          "The given number of control points in each dimension doesn't fit the length of the control point vector.");
    }
    for (auto &&cp : control_points) {
      if (cp.GetDimension() != dimension_) {
        throw std::runtime_error("The dimension has to be the same for all control points.");
      }
      for (int i = 0; i < dimension_; ++i) {
        control_points_.emplace_back(cp.GetValue(i));
      }
    }
  }

  virtual baf::ControlPoint GetControlPoint(std::array<int, DIM> indices) const {
    std::vector<double> coordinates;
    util::MultiIndexHandler<DIM> point_handler = util::MultiIndexHandler<DIM>(number_of_points_);
    point_handler.SetIndices(indices);
    int first = dimension_ * point_handler.Get1DIndex();
    for (int coordinate = 0; coordinate < dimension_; coordinate++) {
      coordinates.push_back(control_points_[first + coordinate]);
    }
    return baf::ControlPoint(coordinates);
  }

  virtual void SetControlPoint(std::array<int, DIM> indices, const baf::ControlPoint &control_point, int dimension) {
    ++number_of_points_[dimension];
    util::MultiIndexHandler<DIM> point_handler = util::MultiIndexHandler<DIM>(number_of_points_);
    point_handler.SetIndices(indices);
    int first = dimension_ * point_handler.Get1DIndex();
    for (int coordinate = 0; coordinate < dimension_; coordinate++) {
      control_points_[first + coordinate] = control_point.GetValue(coordinate);
    }
    --number_of_points_[dimension];
  }

  void AddControlPoints(int number) {
    for (int i = 0; i < number; ++i) {
      for (int j = 0; j < dimension_; ++j) {
        control_points_.emplace_back(0.0);
      }
    }
  }

  int GetNumberOfControlPoints() const {
    return static_cast<int>(control_points_.size()) / dimension_;
  }

  std::array<int, DIM> GetNumberOfPointsInEachDirection() const {
    return number_of_points_;
  }

  std::array<int, DIM> GetMaximumPointIndexInEachDirection() const {
    std::array<int, DIM> maximum_index = number_of_points_;
    for (auto &index : maximum_index) {
      --index;
    }
    return maximum_index;
  }

  void IncrementNumberOfPoints(int dimension) {
    ++number_of_points_[dimension];
  }

  int GetDimension() const {
    return dimension_;
  }

  double GetExpansion() const {
    double expansion = 0;
    for (auto &cp : control_points_) {
      if (cp > expansion) expansion = cp;
    }
    return expansion;
  }

  std::vector<double> GetControlPoints() const {
    return control_points_;
  }

  virtual std::vector<double> GetWeights() const {
    std::vector<double> weights;
    int numberOfWeights = 1;
    for (int i = 0; i < DIM; ++i) {
        numberOfWeights *= number_of_points_[i];
    }
    for (int i = 0; i < numberOfWeights; ++i) {
      weights.emplace_back(1.0);
    }
    return weights;
  }

  virtual double GetWeight(std::array<int, DIM>) const {
    return 1.0;
  }

  std::vector<baf::ControlPoint> GetSplittedControlPoints(int first, int length, int dimension) {
    std::vector<baf::ControlPoint> points;
    std::array<int, DIM> point_handler_length = GetNumberOfPointsInEachDirection();
    point_handler_length[dimension] = length;
    util::MultiIndexHandler<DIM> point_handler(point_handler_length);
    for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
      auto indices = point_handler.GetIndices();
      indices[dimension] += first;
      points.emplace_back(GetControlPoint(indices));
    }
    return points;
  }

 protected:
  int dimension_;
  std::array<int, DIM> number_of_points_;
  std::vector<double> control_points_;
};

}  // namespace spl

#endif  // SRC_SPL_PHYSICAL_SPACE_H_
