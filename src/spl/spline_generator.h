/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

namespace spl {
class SplineGenerator {
 public:
  SplineGenerator() = default;
  virtual ~SplineGenerator() = default;

  virtual ParameterSpace GetParameterSpace() = 0;
  virtual PhysicalSpace GetPhysicalSpace() = 0;
  virtual std::vector<baf::ControlPoint> GetControlPoints() = 0;
};
}

#ifndef SPLINELIB_SPLINE_GENERATOR_H
#define SPLINELIB_SPLINE_GENERATOR_H

#endif //SPLINELIB_SPLINE_GENERATOR_H
