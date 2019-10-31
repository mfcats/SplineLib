/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University
This file is part of SplineLib.
SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.
SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_SQUARE_GENERATOR_H_
#define SRC_SPL_SQUARE_GENERATOR_H_

#include <array>
#include <vector>

#include "src/spl/b_spline.h"
#include "src/baf/knot_vector.h"

namespace splinelib::src::spl {
class SquareGenerator {
 public:
  SquareGenerator();
  SquareGenerator(Degree degree, u_int64_t number_of_knots);

  std::unique_ptr<BSpline<2>> CreateSquare() const;

 private:
  spl::ParameterSpace<2> GenerateParameterSpace() const;

  spl::PhysicalSpace<2> GeneratePhysicalSpace() const;

  Degree degree_;
  u_int64_t number_of_knots_;  // NOLINT
  ParametricCoordinate one_{1};
  ParametricCoordinate zero_{0};
};
}  // namespace splinelib::src::spl

#endif  // SRC_SPL_SQUARE_GENERATOR_H_
