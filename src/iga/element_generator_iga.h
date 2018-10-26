/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_ELEMENT_GENERATOR_IGA_H_
#define SRC_IGA_ELEMENT_GENERATOR_IGA_H_

#include <vector>

#include "element_iga.h"
#include "named_type.h"
#include "spline.h"

namespace iga {
class ElementGenerator {
 public:
  explicit ElementGenerator(std::shared_ptr<spl::Spline<2>> spline);

  std::vector<iga::Element> GetElementList(int i);

  std::vector<int> GetKnotMultiplicity(int i);

  std::vector<ParamCoord> GetInternalKnots(int i);

 private:
  ParamCoord GetLowerElementBound(uint64_t currentKnot, int i);

  ParamCoord GetHigherElementBound(uint64_t currentKnot, int i);

  std::vector<ParamCoord> GetElementNodes(uint64_t currentKnot, int i);

  std::shared_ptr<spl::Spline<2>> spline_;
};
}  // namespace iga

#endif  // SRC_IGA_ELEMENT_GENERATOR_IGA_H_
