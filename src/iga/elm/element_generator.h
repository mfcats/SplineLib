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
#include "multi_index_handler.h"
#include "named_type.h"
#include "spline.h"

namespace iga {
namespace elm {
template<int DIM>
class ElementGenerator {
 public:
  explicit ElementGenerator(std::shared_ptr<spl::Spline<DIM>> spl) : spline_(std::move(spl)) {}

  std::vector<iga::elm::Element> GetElementList(int dir) const {
    std::vector<iga::elm::Element> elements;
    for (uint64_t k = 0; k < spline_->GetKnotVector(dir)->GetNumberOfKnots() - spline_->GetDegree(dir).get() - 1;
         ++k) {
      if ((GetLowerElementBound(k, dir).get() - GetHigherElementBound(k, dir).get()) != 0) {
        elements.emplace_back(Element(1, GetElementNodes(k, dir)));
      }
    }
    return elements;
  }

  std::vector<int> GetKnotMultiplicity(int dir) const {
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

  int GetElementNumberAtParamCoord(std::array<ParamCoord, 2> param_coord) const {
    uint64_t element_number_xi = 0;
    uint64_t element_number_eta = 0;
    std::vector<ParamCoord> unique_knots_xi = GetUniqueKnots(0);
    std::vector<ParamCoord> unique_knots_eta = GetUniqueKnots(1);
    for (uint64_t i = 0; i < unique_knots_xi.size() - 1; ++i) {
      if ((unique_knots_xi[i].get() <= param_coord[1].get()) &&
          (unique_knots_xi[i + 1].get() > param_coord[1].get())) {
        element_number_xi = i;
      }
    }
    for (uint64_t i = 0; i < unique_knots_eta.size() - 1; ++i) {
      if ((unique_knots_eta[i].get() <= param_coord[1].get()) &&
          (unique_knots_eta[i + 1].get() > param_coord[1].get())) {
        element_number_eta = i;
      }
    }
    return Get1DElementIndex(element_number_xi, element_number_eta);
  }

  std::array<int, 2> Get2DElementIndices(int element_number) const {
    element_number += 1;
    int number_of_elements_xi = static_cast<int>(GetElementList(0).size());
    int q = element_number / number_of_elements_xi;
    int r = element_number % number_of_elements_xi;
    std::array<int, 2> element_indices_2d{};
    if (r == 0) {
      element_indices_2d[1] = q - 1;
      element_indices_2d[0] = number_of_elements_xi - 1;
    } else {
      element_indices_2d[1] = q;
      element_indices_2d[0] = r - 1;
    }
    return element_indices_2d;
  }

 private:
  std::vector<ParamCoord> GetInternalKnots(int dir) const {
    std::vector<ParamCoord> knots = spline_->GetKnots()[dir];
    std::vector<ParamCoord> internal_knots;
    int degree = spline_->GetDegree(dir).get();
    for (auto j = uint64_t(degree); j < knots.size() - degree; ++j) {
      internal_knots.emplace_back(knots[j]);
    }
    return internal_knots;
  }

  std::vector<ParamCoord> GetUniqueKnots(int dir) const {
    std::vector<ParamCoord> internal_knots = GetInternalKnots(dir);
    std::vector<ParamCoord> unique_knots;
    for (uint64_t i = 0; i < internal_knots.size() - 1; ++i) {
      if (internal_knots[i].get() != internal_knots[i + 1].get()) {
        unique_knots.emplace_back(internal_knots[i]);
      }
    }
    unique_knots.emplace_back(internal_knots[internal_knots.size() - 1]);
    return unique_knots;
  }

  int Get1DElementIndex(uint64_t element_number_xi, uint64_t element_number_eta) const {
    util::MultiIndexHandler<2> multi_index_handler({static_cast<int>(GetElementList(0).size()),
                                                    static_cast<int>(GetElementList(1).size())});
    multi_index_handler.SetIndices({static_cast<int>(element_number_xi), static_cast<int>(element_number_eta)});
    return multi_index_handler.Get1DIndex();
  }

  ParamCoord GetLowerElementBound(uint64_t currentKnot, int dir) const {
    return spline_->GetKnotVector(dir)->GetKnots()[currentKnot];
  }

  ParamCoord GetHigherElementBound(uint64_t currentKnot, int dir) const {
    return spline_->GetKnotVector(dir)->GetKnots()[currentKnot + 1];
  }

  std::vector<ParamCoord> GetElementNodes(uint64_t currentKnot, int dir) const {
    return {GetLowerElementBound(currentKnot, dir), GetHigherElementBound(currentKnot, dir)};
  }

  std::shared_ptr<spl::Spline<DIM>> spline_;
};
}  // namespace elm
}  // namespace iga

#endif  // SRC_IGA_ELM_ELEMENT_GENERATOR_H_
