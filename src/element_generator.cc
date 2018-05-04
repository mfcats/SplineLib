/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "element_generator.h"

ElementGenerator::ElementGenerator(int degree, const KnotVector &knot_vector)
    : degree_(degree), knot_vector_(knot_vector) {}

std::vector<Element> ElementGenerator::GetElementList() {
  std::vector<Element> elements;
  for (uint64_t currentKnot = 0; currentKnot < knot_vector_.Size() - degree_ - 1; currentKnot++) {
    if ((GetLowerElementBound(currentKnot) - GetHigherElementBound(currentKnot))!=0) {
      elements.emplace_back(Element(1, GetElementNodes(currentKnot)));
    }
  }
  return elements;
}

double ElementGenerator::GetLowerElementBound(uint64_t currentKnot) {
  return knot_vector_[currentKnot];
}

double ElementGenerator::GetHigherElementBound(uint64_t currentKnot) {
  return knot_vector_[currentKnot + 1];
}

std::vector<double> ElementGenerator::GetElementNodes(uint64_t currentKnot) {
  return {GetLowerElementBound(currentKnot), GetHigherElementBound(currentKnot)};
}
