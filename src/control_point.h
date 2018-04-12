/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_CONTROL_POINT_H
#define SPLINELIB_CONTROL_POINT_H

#include <initializer_list>
#include <vector>

#include "definitions.h"

class ControlPoint {
 public:
  explicit ControlPoint(std::initializer_list<double> coordinates);
  explicit ControlPoint(const std::vector<double> &coordinates);

  Dimension dimension() const;
  double GetValue(Dimension dimension) const;

 protected:
  std::vector<double> coordinates_;
};

#endif //SPLINELIB_CONTROL_POINT_H
