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
#include <math.h>

#include "element_generator.h"
#include "src/spl/nurbs.h"

namespace iga {
template<int PARAMETRIC_DIMENSIONALITY>
class MappingHandler {
 public:
  explicit MappingHandler(std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>> spl) : spline_(std::move(spl)) {}

  arma::dmat GetDxiDx(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord) const {
    arma::dmat dx_dxi_sq = GetDxDxi(param_coord).submat(0, 0, static_cast<uint64_t>(PARAMETRIC_DIMENSIONALITY - 1),
                                                        static_cast<uint64_t>(PARAMETRIC_DIMENSIONALITY - 1));
    return dx_dxi_sq.i();
  }

  double GetJacobianDeterminant(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord) const {
    return pow(abs(arma::det(GetDxDxitilde(param_coord).t() * GetDxDxitilde(param_coord))), 0.5);
  }

  std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> Reference2ParameterSpace(int element_number, std::array<double, PARAMETRIC_DIMENSIONALITY> itg_pnts) const {
    iga::elm::ElementGenerator<PARAMETRIC_DIMENSIONALITY> elm_gen(spline_);
    std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coords{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      iga::elm::Element elm = elm_gen.GetElementList(i)[elm_gen.GetElementIndices(element_number)[i]];
      param_coords[i] = ParametricCoordinate{((elm.GetUpperBound() - elm.GetLowerBound()).Get() * itg_pnts[i] +
          (elm.GetUpperBound() + elm.GetLowerBound()).Get()) / 2.0};
    }
    return param_coords;
  }

 private:
  arma::dmat GetDxDxitilde(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord) const {
    return GetDxDxi(param_coord) * GetDxiDxitilde(param_coord);
  }

  arma::dmat GetDxDxi(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord) const {
    int cp_dim = spline_->GetPointDim();
    arma::dmat dx_dxi(static_cast<uint64_t>(cp_dim), static_cast<uint64_t>(PARAMETRIC_DIMENSIONALITY), arma::fill::zeros);
    for (int i = 0; i < cp_dim; ++i) {
      for (int j = 0; j < PARAMETRIC_DIMENSIONALITY; ++j) {
        std::array<int, PARAMETRIC_DIMENSIONALITY> derivative{};
        derivative[j] = 1;
        dx_dxi(static_cast<uint64_t>(i), static_cast<uint64_t>(j)) =
            spline_->EvaluateDerivative(param_coord, {i}, derivative)[0];
      }
    }
    return dx_dxi;
  }

  arma::dmat GetDxiDxitilde(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord) const {
    arma::dmat dxi_dxitilde(static_cast<uint64_t>(PARAMETRIC_DIMENSIONALITY), static_cast<uint64_t>(PARAMETRIC_DIMENSIONALITY), arma::fill::zeros);
    std::array<size_t, PARAMETRIC_DIMENSIONALITY> knot_span{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      knot_span[i] = static_cast<size_t>(spline_->GetKnotVector(i)->GetKnotSpan(param_coord[i]).Get());
      dxi_dxitilde(static_cast<uint64_t>(i), static_cast<uint64_t>(i)) =
          (spline_->GetKnotVector(i)->GetKnot(knot_span[i] + 1).Get() -
              spline_->GetKnotVector(i)->GetKnot(knot_span[i]).Get()) / 2;
    }
    return dxi_dxitilde;
  }

  std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>> spline_;
};
}  // namespace iga

#endif  // SRC_IGA_MAPPING_HANDLER_H_
