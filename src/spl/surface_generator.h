/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_SURFACE_GENERATOR_H_
#define SRC_SPL_SURFACE_GENERATOR_H_

#include <vector>
#include <cmath>

#include "nurbs_generator.h"
#include "nurbs.h"

namespace spl {
class SurfaceGenerator : public NURBSGenerator<2> {
 public:
  SurfaceGenerator(std::shared_ptr<NURBS<1>> const &nurbs_T, std::shared_ptr<NURBS<1>> const &nurbs_C);

  SurfaceGenerator(std::shared_ptr<NURBS<1>> const &nurbs_T,
      std::shared_ptr<NURBS<1>> const &nurbs_C,
      int nbInterpolations, std::vector<double> scaling);

  std::shared_ptr<ParameterSpace<2>> JoinParameterSpaces(std::shared_ptr<ParameterSpace<1>> const &space_1,
      std::shared_ptr<ParameterSpace<1>> const &space_2) const;

  std::shared_ptr<WeightedPhysicalSpace<2>> JoinPhysicalSpaces(std::shared_ptr<PhysicalSpace<1>> const &space_1,
      std::shared_ptr<PhysicalSpace<1>> const &space_2) const;

  double ComputeNormal(
      std::shared_ptr<spl::NURBS<1>> const &nurbs, ParamCoord param_coord, int direction) const;

  std::array<double, 3> CrossProduct( std::vector<double> a, std::vector<double> b) const;

  std::array<double, 3> CrossProduct( std::array<double, 3> a, std::array<double, 3> b) const;

  baf::KnotVector AverageKnots(Degree degree, int nbControlPoints);
};
}  // namespace spl

#endif  // SRC_SPL_SURFACE_GENERATOR_H_
