/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_UTIL_ELEMENT_GENERATOR_H_
#define SRC_UTIL_ELEMENT_GENERATOR_H_

#include <vector>

#include "src/spl/spline.h"
#include "src/util/element.h"
#include "src/util/multi_index_handler.h"

namespace splinelib::src::util {
template<int PARAMETRIC_DIMENSIONALITY>
class ElementGenerator {
 public:
  explicit ElementGenerator(std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>> spl) : spl_(std::move(spl)) {
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      for (uint64_t j = 0; j < spl_->GetKnotVector(i)->GetNumberOfKnots() - spl_->GetDegree(i).Get() - 1; ++j) {
        if ((spl_->GetKnotVector(i)->GetKnot(j).Get() - spl_->GetKnotVector(i)->GetKnot(j + 1).Get()) != 0) {
          elements_[i].emplace_back(Element({spl_->GetKnotVector(i)->GetKnot(j),
                                             spl_->GetKnotVector(i)->GetKnot(j + 1)}));
        }
      }
    }
  }

  std::vector<util::Element> GetElementList(int dir) const {
    return elements_[dir];
  }

 private:
  std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>> spl_;
  std::array<std::vector<util::Element>, PARAMETRIC_DIMENSIONALITY> elements_;
};
}  // namespace splinelib::src::util

#endif  // SRC_UTIL_ELEMENT_GENERATOR_H_
