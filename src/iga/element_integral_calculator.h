/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_ELEMENT_INTEGRAL_CALCULATOR_H_
#define SRC_IGA_ELEMENT_INTEGRAL_CALCULATOR_H_

#include <armadillo>
#include <vector>

#include "basis_function_handler.h"
#include "connectivity_handler.h"
#include "element_integration_point.h"
#include "integration_point.h"
#include "mapping_handler.h"
#include "nurbs.h"

namespace iga {
class ElementIntegralCalculator {
 public:
  explicit ElementIntegralCalculator(std::shared_ptr<spl::NURBS<2>> spl);

  void GetLaplaceElementIntegral(int element_number, const iga::itg::IntegrationRule &rule,
                                 const std::shared_ptr<arma::dmat> &matA) const;

  void GetLaplaceElementIntegral(int element_number, const iga::itg::IntegrationRule &rule,
                                 const std::shared_ptr<arma::dvec> &vecB,
                                 const std::shared_ptr<arma::dvec> &srcCp) const;

 private:
  std::shared_ptr<spl::NURBS<2>> spline_;
  std::shared_ptr<iga::BasisFunctionHandler> baf_handler_;
  std::shared_ptr<iga::ConnectivityHandler> connectivity_handler_;
};
}  // namespace iga

#endif  // SRC_IGA_ELEMENT_INTEGRAL_CALCULATOR_H_
