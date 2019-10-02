/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_CONNECTIVITY_HANDLER_H_
#define SRC_IGA_CONNECTIVITY_HANDLER_H_

#include <vector>

#include "element_generator.h"
#include "multi_index_handler.h"
#include "spline.h"

namespace iga {
template<int PARAMETRIC_DIMENSIONALITY>
class ConnectivityHandler {
 public:
  explicit ConnectivityHandler(std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>> spl) : spline_(std::move(spl)) {
    elm_gen_ = std::make_shared<iga::elm::ElementGenerator<PARAMETRIC_DIMENSIONALITY>>(spline_);
  }

  int GetGlobalIndex(int element_number, int local_index) {
    return Get1DGlobalIndex(GetGlobalIndicesPerParametricDirection(element_number, local_index));
  }

 private:
  std::array<int, PARAMETRIC_DIMENSIONALITY> GetGlobalIndicesPerParametricDirection(int elm_num, int local_index) {
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> mult_ind_handl_elm(elm_gen_->GetNumElementsPerParamDir());
    mult_ind_handl_elm = mult_ind_handl_elm + elm_num;
    std::array<int, PARAMETRIC_DIMENSIONALITY> num_non_zero_baf{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      num_non_zero_baf[i] = spline_->GetDegree(i).get() + 1;
    }
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> mult_ind_handl_baf(num_non_zero_baf);
    mult_ind_handl_baf = mult_ind_handl_baf + local_index;
    std::array<int, PARAMETRIC_DIMENSIONALITY> knot_mult_index_shift = elm_gen_->GetKnotMultiplicityIndexShift(elm_num);
    std::array<int, PARAMETRIC_DIMENSIONALITY> global_indices{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      global_indices[i] = mult_ind_handl_baf[i] + mult_ind_handl_elm[i] + knot_mult_index_shift[i];
    }
    return global_indices;
  }

  int Get1DGlobalIndex(std::array<int, PARAMETRIC_DIMENSIONALITY> global_indices) {
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> mult_ind_handl_cp(spline_->GetPointsPerDirection());
    mult_ind_handl_cp.SetIndices(global_indices);
    return mult_ind_handl_cp.Get1DIndex() + 1;
  }

  std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>> spline_;
  std::shared_ptr<iga::elm::ElementGenerator<PARAMETRIC_DIMENSIONALITY>> elm_gen_;
};
}  // namespace iga

#endif  // SRC_IGA_CONNECTIVITY_HANDLER_H_
