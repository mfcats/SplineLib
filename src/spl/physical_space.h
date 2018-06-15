/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_PHYSICAL_SPACE_H
#define SRC_SPL_PHYSICAL_SPACE_H

#include <stdexcept>
#include <vector>

#include "control_point.h"
#include "multi_index_handler.h"

namespace spl {
template<int DIM>
class PhysicalSpace {
 public:
  PhysicalSpace() = default;
  explicit PhysicalSpace(const std::vector<baf::ControlPoint> &control_points, std::array<int, DIM> number_of_points)
      : dimension_(control_points[0].GetDimension()), point_handler_(number_of_points) {
    int total_number_of_points = 1;
    for (int dim = 0; dim < DIM; dim++) {
      total_number_of_points *= number_of_points[dim];
    }
    if (total_number_of_points != control_points.size()) {
      auto a = control_points.size();
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

  int GetDimension() const {
    return dimension_;
  }

  baf::ControlPoint GetControlPoint(std::array<int, DIM> indices) {
    std::vector<double> coordinates;
    point_handler_.SetIndices(indices);
    int first = dimension_ * point_handler_.Get1DIndex();
    for (int coordinate = 0; coordinate < dimension_; coordinate++) {
      coordinates.push_back(control_points_[first + coordinate]);
    }
    return baf::ControlPoint(coordinates);
  }

  double GetPoint(int point) const {
    return control_points_[point];
  }

  int GetNumberOfControlPoints() const {
    return static_cast<int>(control_points_.size()) / dimension_;
  }

 private:
  std::vector<double> control_points_;
  util::MultiIndexHandler<DIM> point_handler_;
  int dimension_;
};

}  // namespace spl

#endif  // SRC_SPL_PHYSICAL_SPACE_H
