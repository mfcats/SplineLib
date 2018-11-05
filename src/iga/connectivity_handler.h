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
  explicit ConnectivityHandler(std::shared_ptr<spl::Spline<2>> spl);

  std::vector<std::vector<int>> GetConnectivity();

 private:
  void SetConnectivityMatrix();

  void SetElementConnectivity();

  void SetGlobalNodePattern();

  std::shared_ptr<spl::Spline<2>> spline;
  std::vector<std::vector<int>> connectivity;
  std::array<std::vector<std::vector<int>>, 2> element_global;
  std::vector<std::vector<int>> global_pattern;
};
}  // namespace iga

#endif  // SRC_IGA_CONNECTIVITY_HANDLER_H_
