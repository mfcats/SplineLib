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
#include "spline.h"

namespace iga {
namespace elm {
template<int DIM>
class ElementGenerator {
 public:
  explicit ElementGenerator(std::shared_ptr<spl::Spline<DIM>> spl) : spl_(std::move(spl)) {
    for (int i = 0; i < DIM; ++i) {
      for (uint64_t j = 0; j < spl_->GetKnotVector(i)->GetNumberOfKnots() - spl_->GetDegree(i).get() - 1; ++j) {
        if ((spl_->GetKnotVector(i)->GetKnot(j).get() - spl_->GetKnotVector(i)->GetKnot(j + 1).get()) != 0) {
          elements_[i].emplace_back(Element(1, {spl_->GetKnotVector(i)->GetKnot(j),
                                                spl_->GetKnotVector(i)->GetKnot(j + 1)}));
        }
      }
    }
  }

  std::vector<iga::elm::Element> GetElementList(int dir) const {
    return elements_[dir];
  }

  std::array<int, DIM> GetNumberOfElements() const {
    std::array<int, DIM> num_elms{};
    for (int i = 0; i < DIM; ++i) {
      num_elms[i] = elements_[i].size();
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
    std::array<int, DIM> element_indices{};
    for (int i = 0; i < DIM; ++i) {
      baf::KnotVector unique_kv(GetUniqueKnots(i));
      element_indices[i] = unique_kv.GetKnotSpan(param_coords[i]).get();
    }
    return Get1DElementIndex(element_indices);
  }

  std::array<int, DIM> GetElementIndices(int element_index) const {
    util::MultiIndexHandler<DIM> mult_ind_handl_elm(GetNumberOfElements());
    mult_ind_handl_elm = mult_ind_handl_elm + element_index;
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
    internal_knots.erase(unique(internal_knots.begin(), internal_knots.end()), internal_knots.end());
    return internal_knots;
  }

  int Get1DElementIndex(std::array<int, DIM> element_indices) const {
    util::MultiIndexHandler<DIM> mult_ind_handl_elm(GetNumberOfElements());
    mult_ind_handl_elm.SetIndices(element_indices);
    return mult_ind_handl_elm.Get1DIndex();
  }

  std::shared_ptr<spl::Spline<DIM>> spl_;
  std::array<std::vector<iga::elm::Element>, DIM> elements_;
};
}  // namespace elm
}  // namespace iga

#endif  // SRC_IGA_ELM_ELEMENT_GENERATOR_H_
