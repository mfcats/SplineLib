/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_MAPPING_HANDLER_H_
#define SRC_IGA_MAPPING_HANDLER_H_

#include <armadillo>

#include "element_generator.h"
#include "spline.h"

namespace iga {
class MappingHandler {
 public:
  explicit MappingHandler(std::shared_ptr<spl::Spline<2>> spl);

  arma::dmat GetDxiDx(std::array<ParamCoord, 2> param_coord) const;

  double GetJacobianDeterminant(std::array<ParamCoord, 2> param_coord) const;

  std::array<ParamCoord, 2> Reference2ParameterSpace(int element_number, double itg_pnt_xi, double itg_pnt_eta) const;

 private:
  arma::dmat GetDxDxitilde(std::array<ParamCoord, 2> param_coord) const;

  arma::dmat GetDxDxi(std::array<ParamCoord, 2> param_coord) const;

  arma::dmat GetDxiDxitilde(std::array<ParamCoord, 2> param_coord) const;

  std::shared_ptr<spl::Spline<2>> spline_;
};
}  // namespace iga

#endif  // SRC_IGA_MAPPING_HANDLER_H_
