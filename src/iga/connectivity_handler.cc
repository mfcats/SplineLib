/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "connectivity_handler.h"

iga::ConnectivityHandler::ConnectivityHandler(std::shared_ptr<spl::Spline<2>> spl) : spline_(std::move(spl)) {
  SetGlobalNodePattern();
  SetElementConnectivity();
  SetConnectivityMatrix();
}

int iga::ConnectivityHandler::GetGlobalIndex(int element_number, int local_index) {
  return connectivity_[element_number][local_index];
}

void iga::ConnectivityHandler::SetConnectivityMatrix() {
  iga::elm::ElementGenerator element_generator_(spline_);
  for (uint64_t i = 0; i < element_generator_.GetElementList(1).size(); ++i) {
    for (uint64_t j = 0; j < element_generator_.GetElementList(0).size(); ++j) {
      std::vector<int> temp;
      for (uint64_t k = 0; k < element_global_[1][i].size(); ++k) {
        for (uint64_t l = 0; l < element_global_[0][j].size(); ++l) {
          temp.emplace_back(global_pattern_[element_global_[1][i][k]][element_global_[0][j][l]]);
        }
      }
      connectivity_.emplace_back(temp);
    }
  }
}

void iga::ConnectivityHandler::SetElementConnectivity() {
  iga::elm::ElementGenerator element_generator_(spline_);
  std::array<int, 2> knot_multiplicity = {0, 0};
  for (int i = 0; i < 2; ++i) {
    for (uint64_t j = 0; j < element_generator_.GetElementList(i).size(); ++j) {
      knot_multiplicity[i] += element_generator_.GetKnotMultiplicity(i)[j];
      std::vector<int> temp;
      for (int k = 1; k < spline_->GetDegree(i).get() + 2; ++k) {
        temp.emplace_back(k + j + knot_multiplicity[i] - 1);
      }
      element_global_[i].emplace_back(temp);
    }
  }
}

void iga::ConnectivityHandler::SetGlobalNodePattern() {
  int l = 0;
  for (int i = 0; i < spline_->GetPointsPerDirection()[0]; ++i) {
    std::vector<int> temp;
    for (int j = 0; j < spline_->GetPointsPerDirection()[1]; ++j) {
      ++l;
      temp.emplace_back(l);
    }
    global_pattern_.emplace_back(temp);
  }
}
