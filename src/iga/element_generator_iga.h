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
template<int DIM>
class ElementGenerator {
 public:
  explicit ElementGenerator(std::shared_ptr<spl::Spline<DIM>> spl) : spline_(std::move(spl)) {}

  std::vector<iga::Element> GetElementList(int index) {
    std::vector<iga::Element> elements;
    for (uint64_t k = 0; k < spline_->GetKnotVector(index)->GetNumberOfKnots() - spline_->GetDegree(index).get() - 1;
         ++k) {
      if ((GetLowerElementBound(k, index).get() - GetHigherElementBound(k, index).get()) != 0) {
        elements.emplace_back(Element(1, GetElementNodes(k, index)));
      }
    }
    return elements;
  }

  std::vector<int> GetKnotMultiplicity(int index) {
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

  std::vector<ParamCoord> GetInternalKnots(int index) {
    std::vector<ParamCoord> knots = spline_->GetKnots().at(uint64_t(index));
    std::vector<ParamCoord> internal_knots;
    int degree = spline_->GetDegree(index).get();
    for (auto j = uint64_t(degree); j < knots.size() - degree; ++j) {
      internal_knots.emplace_back(knots[j]);
    }
    return internal_knots;
  }

 private:
  ParamCoord GetLowerElementBound(uint64_t currentKnot, int index) {
    return spline_->GetKnotVector(index)->GetKnots()[currentKnot];
  }

  ParamCoord GetHigherElementBound(uint64_t currentKnot, int index) {
    return spline_->GetKnotVector(index)->GetKnots()[currentKnot + 1];
  }

  std::vector<ParamCoord> GetElementNodes(uint64_t currentKnot, int index) {
    return {GetLowerElementBound(currentKnot, index), GetHigherElementBound(currentKnot, index)};
  }

  std::shared_ptr<spl::Spline<DIM>> spline_;
};
}  // namespace iga

#endif  // SRC_IGA_ELEMENT_GENERATOR_IGA_H_
