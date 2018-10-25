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

#include "element_generator_iga.h"
#include "matrix_utils.h"
#include "spline.h"

#include <iostream>

namespace iga {
class ConnectivityHandler {
 public:
  explicit ConnectivityHandler(const std::shared_ptr<spl::Spline<2>> &spline) {

    SetGlobalNodePattern(spline);
    SetElementConnectivity(spline);

    iga::MatrixUtils::PrintMatrix(element_global[0]);

  }

  void SetElementConnectivity(const std::shared_ptr<spl::Spline<2>> &spline) {
    iga::ElementGenerator element_generator_(spline);
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

  void SetGlobalNodePattern(const std::shared_ptr<spl::Spline<2>> &spline) {
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

 private:
  // global_pattern[i][j] is the i-th row and j-th column of the matrix
  std::vector<std::vector<int>> global_pattern;

  // i-th row, j-th column -> global index of j-th control point in i-th element
  std::array<std::vector<std::vector<int>>, 2> element_global;



};
}  // namespace iga

#endif  // SRC_IGA_CONNECTIVITY_HANDLER_H_
