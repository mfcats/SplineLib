/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_SPLINE_GENERATOR_H_
#define SRC_SPL_SPLINE_GENERATOR_H_

#include "parameter_space.h"
#include "physical_space.h"
#include "weighted_physical_space.h"

namespace spl {
template<int DIM>
class SplineGenerator {
 public:
  SplineGenerator() = default;
  virtual ~SplineGenerator() = default;

  SplineGenerator(KnotVectors<DIM> knot_vector, std::array<Degree, DIM> degree) {
    parameter_space_ = std::make_shared<ParameterSpace<DIM>>(ParameterSpace<DIM>(knot_vector, degree));
  }

  std::shared_ptr<ParameterSpace<DIM>> GetParameterSpace() {
    return parameter_space_;
  }

 protected:
  std::shared_ptr<ParameterSpace<DIM>> parameter_space_;
};
}  // namespace spl

#endif  // SRC_SPL_SPLINE_GENERATOR_H_
