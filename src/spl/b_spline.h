/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_B_SPLINE_H_
#define SRC_B_SPLINE_H_

#include <array>
#include <vector>

#include "spline.h"

namespace spl {
template<int DIM>
class BSpline : public Spline<DIM> {
 public:
  BSpline(const std::array<baf::KnotVector, DIM> &knot_vector,
          std::array<int, DIM> degree,
          const std::vector<baf::ControlPoint> &control_points) : Spline<DIM>(knot_vector, degree, control_points) {}

 private:
  std::vector<double> EvaluateAllNonZeroBasisFunctions(std::array<double, DIM> param_coord) const override {
    auto first_non_zero = this->CreateArrayFirstNonZeroBasisFunction(param_coord);
    auto total_length = this->ArrayTotalLength();
    auto M = this->MultiIndexHandlerShort();

    util::MultiIndexHandler<DIM> multiIndexHandler(total_length);

    std::vector<double> vector(M, 1);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < DIM; ++j) {
        vector[i] *= (*(first_non_zero[j] + multiIndexHandler[j]))->Evaluate(param_coord[j]);
      }
      multiIndexHandler++;
    }
    return vector;
  }
};
}

#endif  // SRC_B_SPLINE_H_
