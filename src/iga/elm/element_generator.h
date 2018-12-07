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
  explicit ElementGenerator(std::shared_ptr<spl::Spline<DIM>> spl) : spl_(std::move(spl)) {}

  std::vector<iga::elm::Element> GetElementList(int dir) const {
    std::vector<iga::elm::Element> elements;
    for (uint64_t k = 0; k < spl_->GetKnotVector(dir)->GetNumberOfKnots() - spl_->GetDegree(dir).get() - 1; ++k) {
      if ((GetLowerElementBound(k, dir).get() - GetHigherElementBound(k, dir).get()) != 0) {
        elements.emplace_back(Element(1, GetElementNodes(k, dir)));
      }
    }
    return elements;
  }

  std::array<int, DIM> GetNumberOfElements() const {
    std::array<int, DIM> num_elms{};
    for (int i = 0; i < DIM; ++i) {
      num_elms[i] = GetElementList(i).size();
    }
    return num_elms;
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

  int GetElementNumberAtParamCoord(std::array<ParamCoord, DIM> param_coords) const {
    std::array<std::vector<ParamCoord>, DIM> unique_knots{};
    std::array<int, DIM> element_indices{};
    for (int i = 0; i < DIM; ++i) {
      unique_knots[i] = GetUniqueKnots(i);
      for (uint64_t j = 0; j < unique_knots[i].size() - 1; ++j) {
        if ((unique_knots[i][j].get() <= param_coords[i].get()) &&
            (unique_knots[i][j + 1].get() > param_coords[i].get())) {
          element_indices[i] = static_cast<int>(j);
        }
      }
    }
    return Get1DElementIndex(element_indices);
  }

  std::array<int, DIM> Get2DElementIndices(int elm_num) const {
    util::MultiIndexHandler<DIM> mult_ind_handl_elm(GetNumberOfElements());
    mult_ind_handl_elm = mult_ind_handl_elm + elm_num;
    return mult_ind_handl_elm.GetIndices();
  }

 private:
  std::vector<ParamCoord> GetInternalKnots(int dir) const {
    std::vector<ParamCoord> knots = spl_->GetKnots()[dir];
    auto first = knots.begin() + spl_->GetDegree(dir).get();
    auto last = knots.end() - spl_->GetDegree(dir).get();
    std::vector<ParamCoord> internal_knots_(first, last);
    return internal_knots_;
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

  int Get1DElementIndex(std::array<int, DIM> element_indices) const {
    util::MultiIndexHandler<DIM> mult_ind_handl_elm(GetNumberOfElements());
    mult_ind_handl_elm.SetIndices(element_indices);
    return mult_ind_handl_elm.Get1DIndex();
  }

  ParamCoord GetLowerElementBound(uint64_t currentKnot, int dir) const {
    return spl_->GetKnotVector(dir)->GetKnots()[currentKnot];
  }

  ParamCoord GetHigherElementBound(uint64_t currentKnot, int dir) const {
    return spl_->GetKnotVector(dir)->GetKnots()[currentKnot + 1];
  }

  std::vector<ParamCoord> GetElementNodes(uint64_t currentKnot, int dir) const {
    return {GetLowerElementBound(currentKnot, dir), GetHigherElementBound(currentKnot, dir)};
  }

  std::shared_ptr<spl::Spline<DIM>> spl_;
};
}  // namespace elm
}  // namespace iga

#endif  // SRC_IGA_ELM_ELEMENT_GENERATOR_H_
