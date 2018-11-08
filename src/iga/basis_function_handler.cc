/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "basis_function_handler.h"

iga::BasisFunctionHandler::BasisFunctionHandler(std::shared_ptr<spl::NURBS<2>> spl) : spline_(std::move(spl)) {
  mapping_handler_ = std::make_shared<iga::MappingHandler>(spline_);
  element_generator_ = std::make_shared<iga::elm::ElementGenerator<2>>(spline_);
}

std::vector<iga::elm::ElementIntegrationPoint> iga::BasisFunctionHandler::EvaluateAllElementNonZeroNURBSBasisFunctions(
    int element_number,
    const iga::itg::IntegrationRule &rule) const {
  std::vector<iga::elm::ElementIntegrationPoint> element_integration_points;
  for (auto &itg_pnt_eta : rule.GetIntegrationPoints()) {
    for (auto &itg_pnt_xi : rule.GetIntegrationPoints()) {
      element_integration_points.emplace_back(iga::elm::ElementIntegrationPoint(
          EvaluateAllNonZeroNURBSBasisFunctions(mapping_handler_->Reference2ParameterSpace(
              element_number, itg_pnt_xi.GetCoordinate(), itg_pnt_eta.GetCoordinate())),
              itg_pnt_xi.GetWeight() * itg_pnt_eta.GetWeight()));
    }
  }
  return element_integration_points;
}

std::vector<iga::elm::ElementIntegrationPoint>
iga::BasisFunctionHandler::EvaluateAllElementNonZeroNURBSBasisFunctionDerivatives(int element_number,
    const iga::itg::IntegrationRule &rule) const {
  std::vector<iga::elm::ElementIntegrationPoint> element_integration_points;
  for (auto &itg_p_eta : rule.GetIntegrationPoints()) {
    for (auto &itg_p_xi : rule.GetIntegrationPoints()) {
      element_integration_points.emplace_back(iga::elm::ElementIntegrationPoint(
          EvaluateAllNonZeroNURBSBasisFunctionDerivatives(mapping_handler_->Reference2ParameterSpace(
              element_number, itg_p_xi.GetCoordinate(), itg_p_eta.GetCoordinate())),
              itg_p_xi.GetWeight() * itg_p_eta.GetWeight()));
    }
  }
  return element_integration_points;
}

std::vector<iga::elm::ElementIntegrationPoint>
iga::BasisFunctionHandler::EvaluateAllElementNonZeroNURBSBafDerivativesPhysical(int element_number,
    const iga::itg::IntegrationRule &rule) const {
  std::vector<iga::elm::ElementIntegrationPoint> element_integration_points;
  for (auto &itg_pnt_eta : rule.GetIntegrationPoints()) {
    for (auto &itg_pnt_xi : rule.GetIntegrationPoints()) {
      std::array<ParamCoord, 2> param_coords = mapping_handler_->Reference2ParameterSpace(
          element_number, itg_pnt_xi.GetCoordinate(), itg_pnt_eta.GetCoordinate());
      element_integration_points.emplace_back(iga::elm::ElementIntegrationPoint(
          EvaluateAllNonZeroNURBSBafDerivativesPhysical(param_coords), itg_pnt_xi.GetWeight() * itg_pnt_eta.GetWeight(),
          mapping_handler_->GetJacobianDeterminant(param_coords)));
    }
  }
  return element_integration_points;
}

std::vector<double> iga::BasisFunctionHandler::EvaluateAllNonZeroNURBSBasisFunctions(
    std::array<ParamCoord, 2> param_coord) const {
  std::array<std::vector<double>, 2> basis_functions = std::array<std::vector<double>, 2>({
    spline_->EvaluateAllNonZeroBasisFunctions(0, param_coord[0]),
    spline_->EvaluateAllNonZeroBasisFunctions(1, param_coord[1])});
  std::vector<double> nurbs_basis_functions;
  double sum = 0;
  int l = 0;
  for (uint64_t i = 0; i < basis_functions[1].size(); ++i) {
    for (uint64_t j = 0; j < basis_functions[0].size(); ++j) {
      double temp = basis_functions[0][j] * basis_functions[1][i] * GetWeight(param_coord, l);
      sum += temp;
      nurbs_basis_functions.emplace_back(temp);
      l += 1;
    }
  }
  for (auto &nbaf : nurbs_basis_functions) {
    nbaf = nbaf / sum;
  }
  return nurbs_basis_functions;
}

std::array<std::vector<double>, 2> iga::BasisFunctionHandler::EvaluateAllNonZeroNURBSBasisFunctionDerivatives(
    std::array<ParamCoord, 2> param_coord) const {
  std::array<std::vector<double>, 2> basis_functions = std::array<std::vector<double>, 2>({
    spline_->EvaluateAllNonZeroBasisFunctions(0, param_coord[0]),
    spline_->EvaluateAllNonZeroBasisFunctions(1, param_coord[1])});
  std::array<std::vector<double>, 2> basis_function_derivatives = std::array<std::vector<double>, 2>({
    spline_->EvaluateAllNonZeroBasisFunctionDerivatives(0, param_coord[0], 1),
    spline_->EvaluateAllNonZeroBasisFunctionDerivatives(1, param_coord[1], 1)});
  std::vector<double> nurbs_basis_functions;
  std::array<std::vector<double>, 2> nurbs_basis_function_derivatives;
  double sum_baf = 0;
  double sum_der_xi = 0;
  double sum_der_eta = 0;
  int l = 0;
  for (uint64_t i = 0; i < basis_function_derivatives[1].size(); ++i) {
    for (uint64_t j = 0; j < basis_function_derivatives[0].size(); ++j) {
      double temp = basis_functions[0][j] * basis_functions[1][i] * GetWeight(param_coord, l);
      double temp_der_xi = basis_function_derivatives[0][j] * basis_functions[1][i] * GetWeight(param_coord, l);
      double temp_der_eta = basis_functions[0][j] * basis_function_derivatives[1][i] * GetWeight(param_coord, l);
      sum_baf += temp;
      sum_der_xi += temp_der_xi;
      sum_der_eta += temp_der_eta;
      nurbs_basis_functions.emplace_back(temp);
      nurbs_basis_function_derivatives[0].emplace_back(temp_der_xi);
      nurbs_basis_function_derivatives[1].emplace_back(temp_der_eta);
      l += 1;
    }
  }
  for (uint64_t i = 0; i < basis_functions.size(); ++i) {
    nurbs_basis_functions[i] = nurbs_basis_functions[i] / sum_baf;
    nurbs_basis_function_derivatives[0][i] = (nurbs_basis_function_derivatives[0][i] * sum_baf -
        nurbs_basis_functions[i] * sum_der_xi) / pow(sum_baf, 2);
    nurbs_basis_function_derivatives[1][i] = (nurbs_basis_function_derivatives[1][i] * sum_baf -
        nurbs_basis_functions[i] * sum_der_eta) / pow(sum_baf, 2);
  }
  return nurbs_basis_function_derivatives;
}

std::array<std::vector<double>, 2> iga::BasisFunctionHandler::EvaluateAllNonZeroNURBSBafDerivativesPhysical(
    std::array<ParamCoord, 2> param_coord) const {
  std::array<std::vector<double>, 2> dr_dx;
  std::array<std::vector<double>, 2> dr_dxi = EvaluateAllNonZeroNURBSBasisFunctionDerivatives(param_coord);
  arma::dmat dxi_dx = mapping_handler_->GetDxiDx(param_coord);
  for (uint64_t i = 0; i < dr_dxi[0].size(); ++i) {
    dr_dx[0].emplace_back(dr_dxi[0][i] * dxi_dx(0, 0) + dr_dxi[1][i] * dxi_dx(1, 0));
    dr_dx[1].emplace_back(dr_dxi[0][i] * dxi_dx(0, 1) + dr_dxi[1][i] * dxi_dx(1, 1));
  }
  return dr_dx;
}

double iga::BasisFunctionHandler::GetWeight(std::array<ParamCoord, 2> param_coord, int local_index) const {
  iga::ConnectivityHandler connectivity_handler(spline_);
  return spline_->GetWeights()[connectivity_handler.GetGlobalIndex(element_generator_->GetElementNumberAtParamCoord(
      param_coord), local_index) - 1];
}
