/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_NURBS_H
#define SPLINELIB_NURBS_H

#include <array>
#include <vector>

#include "b_spline.h"
#include "spline.h"

namespace spl {
template<int DIM>
class NURBS : public Spline<DIM> {
 public:
  NURBS(const std::array<baf::KnotVector, DIM> &knot_vector,
        std::array<int, DIM> degree,
        const std::vector<baf::ControlPoint> &control_points, const std::vector<double> weights) : Spline<DIM>(
      knot_vector,
      degree,
      control_points), weights_(weights) {}

 private:
  std::vector<double> EvaluateAllNonZeroBasisFunctions(std::array<double, DIM> param_coord) const override {
    std::vector<baf::ControlPoint> weights;
    for (int control_point = 0; control_point < weights_.size(); ++control_point) {
      weights.emplace_back(baf::ControlPoint({weights_[control_point]}));
    }

    auto first_non_zero = this->CreateArrayFirstNonZeroBasisFunction(param_coord);
    auto total_length = this->ArrayTotalLength();
    auto M = this->MultiIndexHandlerShort();

    util::MultiIndexHandler<DIM> multiIndexHandler(total_length);

    std::vector<double> vector(M, 1);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < DIM; ++j) {
        auto a = (*(first_non_zero[j] + multiIndexHandler[j]))->Evaluate(param_coord[j]);
        std::cout << a << std::endl;
        vector[i] *= (*(first_non_zero[j] + multiIndexHandler[j]))->Evaluate(param_coord[j])
            * weights_[this->GetKnotVector(0).GetKnotSpan(param_coord[0]) - this->GetDegree(0) + multiIndexHandler[j]];
      }
      auto b = vector[i];
      std::cout << b << std::endl;
      vector[i] /= spl::BSpline<DIM>(std::array<baf::KnotVector, DIM>{this->GetKnotVector(0)},
                                     std::array<int, DIM>{this->GetDegree(0)},
                                     weights).Evaluate(param_coord, {0})[0];
      b = vector[i];
      std::cout << b << std::endl;
      multiIndexHandler++;
    }

    return vector;
  }

  std::vector<double> weights_;
};
} //namespace spl

#endif //SPLINELIB_NURBS_H
