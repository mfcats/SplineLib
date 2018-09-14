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

elm::ElementGenerator::ElementGenerator(int degree, baf::KnotVector knot_vector)
    : degree_(degree), knot_vector_(std::move(knot_vector)) {}

std::vector<elm::Element> elm::ElementGenerator::GetElementList() {
  std::vector<Element> elements;
  for (uint64_t currentKnot = 0; currentKnot < knot_vector_.GetNumberOfKnots() - degree_ - 1; currentKnot++) {
    if ((GetLowerElementBound(currentKnot).get() - GetHigherElementBound(currentKnot).get()) != 0) {
      elements.emplace_back(Element(1, GetElementNodes(currentKnot)));
    }
  }
  return elements;
}

ParamCoord elm::ElementGenerator::GetLowerElementBound(uint64_t currentKnot) {
  return knot_vector_[currentKnot];
}

ParamCoord elm::ElementGenerator::GetHigherElementBound(uint64_t currentKnot) {
  return knot_vector_[currentKnot + 1];
}

std::vector<ParamCoord> elm::ElementGenerator::GetElementNodes(uint64_t currentKnot) {
  return {GetLowerElementBound(currentKnot), GetHigherElementBound(currentKnot)};
}
