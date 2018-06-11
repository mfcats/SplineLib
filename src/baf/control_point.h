/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_BAF_CONTROL_POINT_H_
#define SRC_BAF_CONTROL_POINT_H_

#include <initializer_list>
#include <vector>

namespace baf {
class ControlPoint {
 public:
  explicit ControlPoint(std::initializer_list<double> coordinates);

  explicit ControlPoint(const std::vector<double> &coordinates);

  int GetDimension() const;

  double GetValue(int dimension) const;

 protected:
  std::vector<double> coordinates_;
};
}

#endif  // SRC_BAF_CONTROL_POINT_H_
