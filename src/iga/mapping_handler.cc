/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "mapping_handler.h"

iga::MappingHandler::MappingHandler(std::shared_ptr<spl::Spline<2>> spl) : spline_(std::move(spl)) {}

arma::dmat iga::MappingHandler::GetDxiDx(std::array<ParamCoord, 2> param_coord) const {
  return GetDxDxi(param_coord).i();
}

double iga::MappingHandler::GetJacobianDeterminant(std::array<ParamCoord, 2> param_coord) const {
  return arma::det(GetDxDxitilde(param_coord));
}

arma::dmat iga::MappingHandler::GetDxDxitilde(std::array<ParamCoord, 2> param_coord) const {
  arma::dmat dx_dxitilde(2, 2, arma::fill::zeros);
  arma::dmat dx_dxi = GetDxDxi(param_coord);
  arma::dmat dxi_dxitilde = GetDxiDxitilde(param_coord);
  for (uint64_t i = 0; i < 2; ++i) {
    for (uint64_t j = 0; j < 2; ++j) {
      dx_dxitilde(i, j) = dx_dxi(i, j) * dxi_dxitilde(j, j);
    }
  }
  return dx_dxitilde;
}

arma::dmat iga::MappingHandler::GetDxDxi(std::array<ParamCoord, 2> param_coord) const {
  arma::dmat dx_dxi(2, 2, arma::fill::zeros);
  dx_dxi(0, 0) = spline_->EvaluateDerivative(param_coord, {0}, {1, 0})[0];
  dx_dxi(0, 1) = spline_->EvaluateDerivative(param_coord, {0}, {0, 1})[0];
  dx_dxi(1, 0) = spline_->EvaluateDerivative(param_coord, {1}, {1, 0})[0];
  dx_dxi(1, 1) = spline_->EvaluateDerivative(param_coord, {1}, {0, 1})[0];
  return dx_dxi;
}

arma::dmat iga::MappingHandler::GetDxiDxitilde(std::array<ParamCoord, 2> param_coord) const {
  arma::dmat dxi_dxitilde(2, 2, arma::fill::zeros);
  std::array<size_t, 2> knot_span{};
  knot_span[0] = static_cast<size_t>(spline_->GetKnotVector(0)->GetKnotSpan(param_coord[0]).get());
  knot_span[1] = static_cast<size_t>(spline_->GetKnotVector(1)->GetKnotSpan(param_coord[1]).get());
  dxi_dxitilde(0, 0) = (spline_->GetKnotVector(0)->GetKnot(knot_span[0] + 1).get()
      - spline_->GetKnotVector(0)->GetKnot(knot_span[0]).get()) / 2;
  dxi_dxitilde(1, 1) = (spline_->GetKnotVector(1)->GetKnot(knot_span[1] + 1).get()
      - spline_->GetKnotVector(1)->GetKnot(knot_span[1]).get()) / 2;
  return dxi_dxitilde;
}

std::array<ParamCoord, 2> iga::MappingHandler::Reference2ParameterSpace(int element_number, double itg_pnt_xi,
    double itg_pnt_eta) const {
  iga::elm::ElementGenerator<2> elm_gen(spline_);
  iga::elm::Element element_xi = elm_gen.GetElementList(0)[elm_gen.GetElementIndices(element_number)[0]];
  iga::elm::Element element_eta = elm_gen.GetElementList(1)[elm_gen.GetElementIndices(element_number)[1]];
  ParamCoord upper_xi = element_xi.GetNode(1);
  ParamCoord lower_xi = element_xi.GetNode(0);
  ParamCoord upper_eta = element_eta.GetNode(1);
  ParamCoord lower_eta = element_eta.GetNode(0);
  return {ParamCoord{((upper_xi - lower_xi).get() * itg_pnt_xi + (upper_xi + lower_xi).get()) / 2.0},
          ParamCoord{((upper_eta - lower_eta).get() * itg_pnt_eta + (upper_eta + lower_eta).get()) / 2.0}};
}
