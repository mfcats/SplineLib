/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_ELEMENT_GENERATOR_H_
#define SRC_ELEMENT_GENERATOR_H_

#include <vector>

#include "element.h"
#include "knot_vector.h"

class ElementGenerator {
 public:
  ElementGenerator(int degree, const KnotVector &knot_vector);

  std::vector<Element> GetElementList();

 private:
  double GetLowerElementBound(uint64_t currentKnot);
  double GetHigherElementBound(uint64_t currentKnot);
  std::vector<double> GetElementNodes(uint64_t currentKnot);

  int degree_;
  KnotVector knot_vector_;
};

#endif  // SRC_ELEMENT_GENERATOR_H_
