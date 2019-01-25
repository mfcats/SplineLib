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
          elements_[i].emplace_back(Element({spl_->GetKnotVector(i)->GetKnot(j),
                                             spl_->GetKnotVector(i)->GetKnot(j + 1)}));
        }
      }
    }
  }

  std::vector<iga::elm::Element> GetElementList(int dir) const {
    return elements_[dir];
  }

  std::array<int, DIM> GetNumElementsPerParamDir() const {
    std::array<int, DIM> num_elms{};
    for (int i = 0; i < DIM; ++i) {
      num_elms[i] = elements_[i].size();
    }
    return num_elms;
  }

  int GetNumberOfElements() const {
    int num_elms = 1;
    for (int i = 0; i < DIM; ++i) {
      num_elms *= elements_[i].size();
    }
    return num_elms;
  }

  std::array<int, DIM> GetKnotMultiplicityIndexShift(int elm_num) const {
    std::array<int, DIM> index_shift{};
    std::array<std::vector<ParamCoord>, DIM> internal = GetInternalKnots();
    std::array<std::vector<ParamCoord>, DIM> unique = GetUniqueKnots();
    for (int i = 0; i < DIM; ++i) {
      ParamCoord lower_bound = unique[i][GetElementIndices(elm_num)[i]];
      internal[i].erase(
          std::find(internal[i].rbegin(), internal[i].rend(), lower_bound).base(), internal[i].rbegin().base());
      unique[i].erase(std::find(unique[i].rbegin(), unique[i].rend(), lower_bound).base(), unique[i].rbegin().base());
      index_shift[i] = internal[i].size() - unique[i].size();
    }
    return index_shift;
  }

  int GetElementNumberAtParamCoord(std::array<ParamCoord, DIM> param_coords) const {
    std::array<int, DIM> element_indices{};
    for (int i = 0; i < DIM; ++i) {
      baf::KnotVector unique_kv(GetUniqueKnots()[i]);
      element_indices[i] = unique_kv.GetKnotSpan(param_coords[i]).get();
    }
    return Get1DElementIndex(element_indices);
  }

  std::array<int, DIM> GetElementIndices(int element_index) const {
    util::MultiIndexHandler<DIM> mult_ind_handl_elm(GetNumElementsPerParamDir());
    mult_ind_handl_elm = mult_ind_handl_elm + element_index;
    return mult_ind_handl_elm.GetIndices();
  }

  int Get1DElementIndex(std::array<int, DIM> element_indices) const {
    util::MultiIndexHandler<DIM> mult_ind_handl_elm(GetNumElementsPerParamDir());
    mult_ind_handl_elm.SetIndices(element_indices);
    return mult_ind_handl_elm.Get1DIndex();
  }

 private:
  std::array<std::vector<ParamCoord>, DIM> GetInternalKnots() const {
    std::array<std::vector<ParamCoord>, DIM> internal_knots;
    for (int i = 0; i < DIM; ++i) {
      std::vector<ParamCoord> knots = spl_->GetKnots()[i];
      auto first = knots.begin() + spl_->GetDegree(i).get();
      auto last = knots.end() - spl_->GetDegree(i).get();
      internal_knots[i] = std::vector<ParamCoord>(first, last);
    }
    return internal_knots;
  }

  std::array<std::vector<ParamCoord>, DIM> GetUniqueKnots() const {
    std::array<std::vector<ParamCoord>, DIM> internal_knots = GetInternalKnots();
    for (int i = 0; i < DIM; ++i) {
      internal_knots[i].erase(unique(internal_knots[i].begin(), internal_knots[i].end()), internal_knots[i].end());
    }
    return internal_knots;
  }

  std::shared_ptr<spl::Spline<DIM>> spl_;
  std::array<std::vector<iga::elm::Element>, DIM> elements_;
};
}  // namespace elm
}  // namespace iga

#endif  // SRC_IGA_ELM_ELEMENT_GENERATOR_H_
