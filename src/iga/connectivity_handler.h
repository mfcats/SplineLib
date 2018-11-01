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
#include "spline.h"

namespace iga {
class ConnectivityHandler {
 public:
  explicit ConnectivityHandler(std::shared_ptr<spl::Spline<2>> spl) : spline(std::move(spl)) {
    SetGlobalNodePattern();
    SetElementConnectivity();
    SetConnectivityMatrix();
  }

  std::vector<std::vector<int>> GetConnectivity() {
    return connectivity;
  }

 private:
  void SetConnectivityMatrix() {
    iga::elm::ElementGenerator element_generator_(spline);
    for (uint64_t i = 0; i < element_generator_.GetElementList(1).size(); ++i) {
      for (uint64_t j = 0; j < element_generator_.GetElementList(0).size(); ++j) {
        std::vector<int> temp;
        for (uint64_t k = 0; k < element_global[1][i].size(); ++k) {
          for (uint64_t l = 0; l < element_global[0][j].size(); ++l) {
            temp.emplace_back(global_pattern[element_global[1][i][k]][element_global[0][j][l]]);
          }
        }
        connectivity.emplace_back(temp);
      }
    }
  }

  void SetElementConnectivity() {
    iga::elm::ElementGenerator element_generator_(spline);
    std::array<int, 2> knot_multiplicity = {0, 0};
    for (int i = 0; i < 2; ++i) {
      for (uint64_t j = 0; j < element_generator_.GetElementList(i).size(); ++j) {
        knot_multiplicity[i] += element_generator_.GetKnotMultiplicity(i)[j];
        std::vector<int> temp;
        for (int k = 1; k < spline->GetDegree(i).get() + 2; ++k) {
          temp.emplace_back(k + j + knot_multiplicity[i] - 1);
        }
        element_global[i].emplace_back(temp);
      }
    }
  }

  void SetGlobalNodePattern() {
    int l = 0;
    for (int i = 0; i < spline->GetPointsPerDirection()[0]; ++i) {
      std::vector<int> temp;
      for (int j = 0; j < spline->GetPointsPerDirection()[1]; ++j) {
        ++l;
        temp.emplace_back(l);
      }
      global_pattern.emplace_back(temp);
    }
  }

  std::shared_ptr<spl::Spline<2>> spline;
  std::vector<std::vector<int>> connectivity;
  std::array<std::vector<std::vector<int>>, 2> element_global;
  std::vector<std::vector<int>> global_pattern;
};
}  // namespace iga

#endif  // SRC_IGA_CONNECTIVITY_HANDLER_H_
