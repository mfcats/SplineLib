/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_UTIL_ELEMENT_H_
#define SRC_UTIL_ELEMENT_H_

#include <memory>
#include <vector>

#include "control_point.h"
#include "knot_vector.h"

namespace util {
class Element {
 public:
  explicit Element(const std::array<ParamCoord, 2> &nodes) : nodes_(nodes) {}

  ParamCoord GetLowerBound() const {
    return nodes_[0];
  }

  ParamCoord GetUpperBound() const {
    return nodes_[1];
  }

 private:
  std::array<ParamCoord, 2> nodes_;
};
}  // namespace util

#endif  // SRC_UTIL_ELEMENT_H_
