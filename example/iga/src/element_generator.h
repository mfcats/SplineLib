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

  std::vector<iga::elm::Element> GetElementList(int dir) const {
    return elements_[dir];
  }

  std::array<int, PARAMETRIC_DIMENSIONALITY> GetNumElementsPerParamDir() const {
    std::array<int, PARAMETRIC_DIMENSIONALITY> num_elms{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      num_elms[i] = elements_[i].size();
    }
    return num_elms;
  }

  int GetNumberOfElements() const {
    int num_elms = 1;
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      num_elms *= elements_[i].size();
    }
    return num_elms;
  }

  std::array<int, PARAMETRIC_DIMENSIONALITY> GetKnotMultiplicityIndexShift(int elm_num) const {
    std::array<int, PARAMETRIC_DIMENSIONALITY> index_shift{};
    std::array<std::vector<ParametricCoordinate>, PARAMETRIC_DIMENSIONALITY> internal = GetInternalKnots();
    std::array<std::vector<ParametricCoordinate>, PARAMETRIC_DIMENSIONALITY> unique = GetUniqueKnots();
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      ParametricCoordinate lower_bound = unique[i][GetElementIndices(elm_num)[i]];
      internal[i].erase(
          std::find(internal[i].rbegin(), internal[i].rend(), lower_bound).base(), internal[i].rbegin().base());
      unique[i].erase(std::find(unique[i].rbegin(), unique[i].rend(), lower_bound).base(), unique[i].rbegin().base());
      index_shift[i] = internal[i].size() - unique[i].size();
    }
    return index_shift;
  }

  int GetElementNumberAtParametricCoordinate(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coords) const {
    std::array<int, PARAMETRIC_DIMENSIONALITY> element_indices{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      baf::KnotVector unique_kv(GetUniqueKnots()[i]);
      element_indices[i] = unique_kv.GetKnotSpan(param_coords[i]).Get();
    }
    return Get1DElementIndex(element_indices);
  }

  std::array<int, PARAMETRIC_DIMENSIONALITY> GetElementIndices(int element_index) const {
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> mult_ind_handl_elm(GetNumElementsPerParamDir());
    mult_ind_handl_elm = mult_ind_handl_elm + element_index;
    return mult_ind_handl_elm.GetCurrentIndex();
  }

  int Get1DElementIndex(std::array<int, PARAMETRIC_DIMENSIONALITY> element_indices) const {
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> mult_ind_handl_elm(GetNumElementsPerParamDir());
    mult_ind_handl_elm.SetCurrentIndex(element_indices);
    return mult_ind_handl_elm.Get1DIndex();
  }

 private:
  std::array<std::vector<ParametricCoordinate>, PARAMETRIC_DIMENSIONALITY> GetInternalKnots() const {
    std::array<std::vector<ParametricCoordinate>, PARAMETRIC_DIMENSIONALITY> internal_knots;
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
//      std::vector<ParametricCoordinate> knots = spl_->GetKnotVector(i)->  //spl_->GetKnots()[i];
      auto first = spl_->GetKnotVector(i)->begin() + spl_->GetDegree(i).Get();
      auto last = spl_->GetKnotVector(i)->end() - spl_->GetDegree(i).Get();
      internal_knots[i] = std::vector<ParametricCoordinate>(first, last);
    }
    return internal_knots;
  }

  std::array<std::vector<ParametricCoordinate>, PARAMETRIC_DIMENSIONALITY> GetUniqueKnots() const {
    std::array<std::vector<ParametricCoordinate>, PARAMETRIC_DIMENSIONALITY> internal_knots = GetInternalKnots();
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      internal_knots[i].erase(unique(internal_knots[i].begin(), internal_knots[i].end()), internal_knots[i].end());
    }
    return internal_knots;
  }

  std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>> spl_;
  std::array<std::vector<iga::elm::Element>, PARAMETRIC_DIMENSIONALITY> elements_;
};
}  // namespace elm
}  // namespace iga

#endif  // SRC_IGA_ELM_ELEMENT_GENERATOR_H_
