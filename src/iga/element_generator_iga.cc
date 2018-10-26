/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "element_generator_iga.h"
#include "spline.h"

iga::ElementGenerator::ElementGenerator(std::shared_ptr<spl::Spline<2>> spl) : spline_(std::move(spl)) {}

std::vector<iga::Element> iga::ElementGenerator::GetElementList(int index) {
  std::vector<iga::Element> elements;
    for (uint64_t k = 0; k < spline_->GetKnotVector(index)->GetNumberOfKnots() - spline_->GetDegree(index).get() - 1;
         ++k) {
      if ((GetLowerElementBound(k, index).get() - GetHigherElementBound(k, index).get()) != 0) {
        elements.emplace_back(Element(1, GetElementNodes(k, index)));
      }
    }
  return elements;
}

std::vector<int> iga::ElementGenerator::GetKnotMultiplicity(int index) {
  std::vector<int> knot_multiplicity;
  std::vector<ParamCoord> internal_knots = GetInternalKnots(index);
  int temp = 0;
  for (uint64_t j = 0; j < internal_knots.size() - 1; ++j) {
    if (internal_knots[j].get() == internal_knots[j + 1].get()) {
      temp += 1;
    } else if ((internal_knots[j].get() != internal_knots[j + 1].get()) || (j == internal_knots.size() - 2)) {
      knot_multiplicity.emplace_back(temp);
      temp = 0;
    }
  }
  return knot_multiplicity;
}

std::vector<ParamCoord > iga::ElementGenerator::GetInternalKnots(int index) {
  std::vector<ParamCoord> knots = spline_->GetKnots().at(uint64_t(index));
  std::vector<ParamCoord> internal_knots;
  int degree = spline_->GetDegree(index).get();
  for (auto j = uint64_t(degree); j < knots.size() - degree; ++j) {
    internal_knots.emplace_back(knots[j]);
  }
  return internal_knots;
}

ParamCoord iga::ElementGenerator::GetLowerElementBound(uint64_t currentKnot, int index) {
  return spline_->GetKnotVector(index)->GetKnots()[currentKnot];
}

ParamCoord iga::ElementGenerator::GetHigherElementBound(uint64_t currentKnot, int index) {
  return spline_->GetKnotVector(index)->GetKnots()[currentKnot + 1];
}

std::vector<ParamCoord> iga::ElementGenerator::GetElementNodes(uint64_t currentKnot, int index) {
  return {GetLowerElementBound(currentKnot, index), GetHigherElementBound(currentKnot, index)};
}
