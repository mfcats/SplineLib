/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_ELM_ELEMENT_GENERATOR_H_
#define SRC_IGA_ELM_ELEMENT_GENERATOR_H_

#include <vector>

#include "element.h"
#include "named_type.h"
#include "spline.h"

namespace iga {
namespace elm {
template<int DIM>
class ElementGenerator {
 public:
  explicit ElementGenerator(std::shared_ptr<spl::Spline<DIM>> spl) : spline_(std::move(spl)) {}

  std::vector<iga::elm::Element> GetElementList(int dir) {
    std::vector<iga::elm::Element> elements;
    for (uint64_t k = 0; k < spline_->GetKnotVector(dir)->GetNumberOfKnots() - spline_->GetDegree(dir).get() - 1;
         ++k) {
      if ((GetLowerElementBound(k, dir).get() - GetHigherElementBound(k, dir).get()) != 0) {
        elements.emplace_back(Element(1, GetElementNodes(k, dir)));
      }
    }
    return elements;
  }

  std::vector<int> GetKnotMultiplicity(int dir) {
    std::vector<int> knot_multiplicity;
    std::vector<ParamCoord> internal_knots = GetInternalKnots(dir);
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

  std::vector<ParamCoord> GetInternalKnots(int dir) {
    std::vector<ParamCoord> knots = spline_->GetKnots().at(uint64_t(dir));
    std::vector<ParamCoord> internal_knots;
    int degree = spline_->GetDegree(dir).get();
    for (auto j = uint64_t(degree); j < knots.size() - degree; ++j) {
      internal_knots.emplace_back(knots[j]);
    }
    return internal_knots;
  }

 private:
  ParamCoord GetLowerElementBound(uint64_t currentKnot, int dir) {
    return spline_->GetKnotVector(dir)->GetKnots()[currentKnot];
  }

  ParamCoord GetHigherElementBound(uint64_t currentKnot, int dir) {
    return spline_->GetKnotVector(dir)->GetKnots()[currentKnot + 1];
  }

  std::vector<ParamCoord> GetElementNodes(uint64_t currentKnot, int dir) {
    return {GetLowerElementBound(currentKnot, dir), GetHigherElementBound(currentKnot, dir)};
  }

  std::shared_ptr<spl::Spline<DIM>> spline_;
};
}  // namespace elm
}  // namespace iga

#endif  // SRC_IGA_ELM_ELEMENT_GENERATOR_H_
