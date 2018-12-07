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
template<int DIM>
class ConnectivityHandler {
 public:
  explicit ConnectivityHandler(std::shared_ptr<spl::Spline<DIM>> spl) : spline_(std::move(spl)) {
    elm_gen_ = std::make_shared<iga::elm::ElementGenerator<2>>(spline_);
    SetConnectivityForParametricDirections();
    SetConnectivityMatrix();
  }

  int GetGlobalIndex(int element_number, int local_index) {
    return connectivity_[element_number][local_index];
  }

 private:
  void SetConnectivityMatrix() {
    for (uint64_t i = 0; i < elm_gen_->GetElementList(1).size(); ++i) {
      for (uint64_t j = 0; j < elm_gen_->GetElementList(0).size(); ++j) {
        std::vector<int> temp;
        for (uint64_t k = 0; k < global_param_[1][i].size(); ++k) {
          for (uint64_t l = 0; l < global_param_[0][j].size(); ++l) {
            temp.emplace_back(Get1DGlobalIndex({global_param_[0][j][l], global_param_[1][i][k]}));
          }
        }
        connectivity_.emplace_back(temp);
      }
    }
  }

  void SetConnectivityForParametricDirections() {
    std::array<int, DIM> knot_multiplicity_sum{};
    for (int i = 0; i < DIM; ++i) {
      for (uint64_t j = 0; j < elm_gen_->GetElementList(i).size(); ++j) {
        knot_multiplicity_sum[i] += elm_gen_->GetKnotMultiplicity(i)[j];
        std::vector<int> temp;
        for (int k = 1; k < spline_->GetDegree(i).get() + 2; ++k) {
          temp.emplace_back(k + j + knot_multiplicity_sum[i] - 1);
        }
        global_param_[i].emplace_back(temp);
      }
    }
  }

  int Get1DGlobalIndex(std::array<int, DIM> global_indices) {
    util::MultiIndexHandler<DIM> mult_ind_handl(spline_->GetPointsPerDirection());
    mult_ind_handl.SetIndices(global_indices);
    return mult_ind_handl.Get1DIndex() + 1;
  }

  std::shared_ptr<spl::Spline<DIM>> spline_;
  std::shared_ptr<iga::elm::ElementGenerator<DIM>> elm_gen_;
  std::vector<std::vector<int>> connectivity_;
  std::array<std::vector<std::vector<int>>, DIM> global_param_;
};
}  // namespace iga

#endif  // SRC_IGA_CONNECTIVITY_HANDLER_H_
