/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_LINEAR_EQUATION_ASSEMBLER_H_
#define SRC_IGA_LINEAR_EQUATION_ASSEMBLER_H_

#include <armadillo>

#include "element_integral_calculator.h"
#include "element_generator.h"
#include "integration_rule.h"
#include "matrix.h"
#include "nurbs.h"

namespace iga {
class LinearEquationAssembler {
 public:
  explicit LinearEquationAssembler(std::shared_ptr<spl::NURBS<2>> spl) : spline_(std::move(spl)) {
    elm_gen_ = std::make_shared<iga::elm::ElementGenerator<2>>(spline_);
  }

  void GetLeftSide(const iga::itg::IntegrationRule &rule, const std::shared_ptr<arma::dmat> &matA,
                   const iga::ElementIntegralCalculator &elm_itg_calc) const {
    int num_elements = static_cast<int>(elm_gen_->GetElementList(0).size() * elm_gen_->GetElementList(1).size());
    for (int e = 0; e < num_elements; ++e) {
      elm_itg_calc.GetLaplaceElementIntegral(e, rule, matA);
    }
  }

 private:
  std::shared_ptr<spl::NURBS<2>> spline_;
  std::shared_ptr<iga::elm::ElementGenerator<2>> elm_gen_;
};
}  // namespace iga

#endif  // SRC_IGA_LINEAR_EQUATION_ASSEMBLER_H_
