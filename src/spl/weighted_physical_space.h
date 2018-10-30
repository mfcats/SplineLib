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

  virtual baf::ControlPoint GetHomogenousControlPoint(std::array<int, DIM> indices) const {
    std::vector<double> coordinates;
    util::MultiIndexHandler<DIM> point_handler = util::MultiIndexHandler<DIM>(this->number_of_points_);
    point_handler.SetIndices(indices);
    int first = this->dimension_ * point_handler.Get1DIndex();
    for (int coordinate = 0; coordinate < this->dimension_; coordinate++) {
      coordinates.push_back(this->control_points_[first + coordinate] * weights_[first / this->dimension_]);
    }
    return baf::ControlPoint(coordinates);
  }

  virtual double GetWeight(std::array<int, DIM> indices) const override {
    util::MultiIndexHandler<DIM> point_handler = util::MultiIndexHandler<DIM>(this->number_of_points_);
    point_handler.SetIndices(indices);
    int first = this->dimension_ * point_handler.Get1DIndex();
    return weights_[first / this->dimension_];
  }

  std::vector<double> GetWeights() const override {
    return weights_;
  }

 private:
  std::vector<double> weights_;
};
}  // namespace spl
#endif  // SRC_SPL_WEIGHTED_PHYSICAL_SPACE_H_
