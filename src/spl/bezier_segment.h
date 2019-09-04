/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_BEZIER_SEGMENT_H_
#define SRC_SPL_BEZIER_SEGMENT_H_

#include <array>
#include <functional>
#include <numeric>
#include <vector>

#include "control_point.h"
#include "knot_vector.h"
#include "multi_index_handler.h"

namespace spl {
template<int DIM>
class BezierSegment {
 public:
  BezierSegment(std::array<Degree, DIM> degrees, const std::vector<baf::ControlPoint> &control_points,
      std::array<bool, DIM> is_bezier_in_direction, std::array<int, DIM> num_control_points)
        : control_points_(control_points), num_control_points_(num_control_points),
          is_bezier_in_direction_(is_bezier_in_direction), degrees_(degrees) {
    uint64_t total_number_of_points = 1;
    for (int i = 0; i < DIM; ++i) {
      total_number_of_points *= num_control_points_[i];
    }
    if (total_number_of_points != control_points.size()) {
      throw std::runtime_error(
          "The given number of control points in each dimension doesn't fit the length of the control point vector.");
    }
  }

  double ReduceDegree(int dimension) {
    if (!is_bezier_in_direction_[dimension]) {
      throw std::runtime_error(
          "The degree in this direction cannot be reduced! The spline is not a Bézier spline in that direction!");
    }
    double max_error_bound = 0.0;
    util::MultiIndexHandler<DIM> cp_handler_old(num_control_points_);
    std::array<Degree, DIM> delta_degrees{};
    delta_degrees[dimension] = Degree{-1};
    SetDegreeAndNumberOfControlPoints(delta_degrees);
    util::MultiIndexHandler<DIM> cp_handler_new(num_control_points_);
    std::vector<baf::ControlPoint> control_points_new(cp_handler_new.Get1DLength(), baf::ControlPoint({0.0}));
    auto indices_new = cp_handler_new.Get1DIndicesForFixedDimension(dimension, 0);
    auto indices_old = cp_handler_old.Get1DIndicesForFixedDimension(dimension, 0);
    for (uint64_t i = 0; i < indices_new.size(); ++i) {
      control_points_new[indices_new[i]] = control_points_[indices_old[i]];
    }
    indices_new = cp_handler_new.Get1DIndicesForFixedDimension(dimension, degrees_[dimension].get());
    indices_old = cp_handler_old.Get1DIndicesForFixedDimension(dimension, degrees_[dimension].get() + 1);
    for (uint64_t i = 0; i < indices_new.size(); ++i) {
      control_points_new[indices_new[i]] = control_points_[indices_old[i]];
    }
    int r = (degrees_[dimension].get()) / 2;
    int end_of_first_loop = r;
    if ((degrees_[dimension].get() + 1) % 2 != 0) end_of_first_loop = r - 1;
    for (int i = 1; i <= end_of_first_loop; ++i) {
      auto indices_new_minus_one = cp_handler_new.Get1DIndicesForFixedDimension(dimension, i - 1);
      indices_new = cp_handler_new.Get1DIndicesForFixedDimension(dimension, i);
      indices_old = cp_handler_old.Get1DIndicesForFixedDimension(dimension, i);
      double alpha = static_cast<double>(i) / static_cast<double>(degrees_[dimension].get() + 1);
      for (uint64_t j = 0; j < indices_new.size(); ++j) {
        control_points_new[indices_new[j]] = (control_points_[indices_old[j]]
            - control_points_new[indices_new_minus_one[j]] * alpha) * (1 / (1 - alpha));
      }
    }
    int lower = static_cast<int>(r + 1);
    int upper = degrees_[dimension].get() - 1;
    for (int i = upper; i >= lower; --i) {
      auto indices_new_plus_one = cp_handler_new.Get1DIndicesForFixedDimension(dimension, i + 1);
      indices_new = cp_handler_new.Get1DIndicesForFixedDimension(dimension, i);
      auto indices_old_plus_one = cp_handler_old.Get1DIndicesForFixedDimension(dimension, i + 1);
      double alpha = static_cast<double>(i + 1) / static_cast<double>(degrees_[dimension].get() + 1);
      for (uint64_t j = 0; j < indices_new.size(); ++j) {
        control_points_new[indices_new[j]] = (control_points_[indices_old_plus_one[j]]
            - control_points_new[indices_new_plus_one[j]] * (1 - alpha)) * (1 / alpha);
      }
    }
    if ((degrees_[dimension].get() + 1) % 2 != 0) {
      auto indices_new_r = cp_handler_new.Get1DIndicesForFixedDimension(dimension, r);
      auto indices_new_r_minus_one = cp_handler_new.Get1DIndicesForFixedDimension(dimension, r - 1);
      auto indices_new_r_plus_one = cp_handler_new.Get1DIndicesForFixedDimension(dimension, r + 1);
      auto indices_old_r = cp_handler_old.Get1DIndicesForFixedDimension(dimension, r);
      auto indices_old_r_plus_one = cp_handler_old.Get1DIndicesForFixedDimension(dimension, r + 1);
      double alpha_r = static_cast<double>(r) / static_cast<double>((degrees_[dimension].get() + 1));
      double alpha_r_plus_one = static_cast<double>(r + 1) / static_cast<double>((degrees_[dimension].get() + 1));
      for (uint64_t j = 0; j < indices_new.size(); ++j) {
        baf::ControlPoint P_r_L = (control_points_[indices_old_r[j]]
            - control_points_new[indices_new_r_minus_one[j]] * alpha_r) * (1 / (1 - alpha_r));
        baf::ControlPoint P_r_R = (control_points_[indices_old_r_plus_one[j]]
            - control_points_new[indices_new_r_plus_one[j]] * (1 - alpha_r_plus_one)) * (1 / alpha_r_plus_one);
        control_points_new[indices_new_r[j]] = (P_r_L + P_r_R) * 0.5;
        if ((P_r_L - P_r_R).GetEuclideanNorm() > max_error_bound) max_error_bound = (P_r_L - P_r_R).GetEuclideanNorm();
      }
    }
    if ((degrees_[dimension].get() + 1) % 2 == 0) {
      auto indices_new_r = cp_handler_new.Get1DIndicesForFixedDimension(dimension, r);
      auto indices_new_r_plus_one = cp_handler_new.Get1DIndicesForFixedDimension(dimension, r + 1);
      auto indices_old_r_plus_one = cp_handler_old.Get1DIndicesForFixedDimension(dimension, r + 1);
      for (uint64_t j = 0; j < indices_new.size(); ++j) {
        double current_max_error_bound =
            (control_points_[indices_old_r_plus_one[j]] - ((control_points_new[indices_new_r[j]]
                + control_points_new[indices_new_r_plus_one[j]]) * 0.5)).GetEuclideanNorm();
        if (current_max_error_bound > max_error_bound) max_error_bound = current_max_error_bound;
      }
    }
    control_points_ = control_points_new;
    return max_error_bound;
  }

  void SetDegreeAndNumberOfControlPoints(std::array<Degree, DIM> delta_degrees) {
    for (int i = 0; i < DIM; ++i) {
      degrees_[i] = degrees_[i] + delta_degrees[i];
      num_control_points_[i] += delta_degrees[i].get();
    }
  }

  baf::ControlPoint GetControlPoint(int index) {
    return control_points_[index];
  }

  Degree GetDegree(int dimension) {
    return degrees_[dimension];
  }

 private:
  std::vector<baf::ControlPoint> control_points_;
  std::array<int, DIM> num_control_points_;
  std::array<bool, DIM> is_bezier_in_direction_;
  std::array<Degree, DIM> degrees_;
};
}  // namespace spl

#endif  // SRC_SPL_BEZIER_SEGMENT_H_
