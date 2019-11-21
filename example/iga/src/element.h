/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_IGA_ELM_ELEMENT_H_
#define SRC_IGA_ELM_ELEMENT_H_

#include <memory>
#include <vector>

#include "src/spl/control_point.h"
#include "src/baf/knot_vector.h"

namespace iga {
namespace elm {
class Element {
 public:
  explicit Element(const std::array<ParametricCoordinate, 2> &nodes) : nodes_(nodes) {}

  ParametricCoordinate GetLowerBound() const {
    return nodes_[0];
  }

  ParametricCoordinate GetUpperBound() const {
    return nodes_[1];
  }

 private:
  std::array<ParametricCoordinate, 2> nodes_;
};
}  // namespace elm
}  // namespace iga

#endif  // SRC_IGA_ELM_ELEMENT_H_
