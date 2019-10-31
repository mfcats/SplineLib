/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_BASIS_FUNCTION_HANDLER_H_
#define SRC_IGA_BASIS_FUNCTION_HANDLER_H_

#include <math.h>
#include <armadillo>
#include <array>
#include <vector>

#include "connectivity_handler.h"
#include "element.h"
#include "element_generator.h"
#include "element_integration_point.h"
#include "integration_rule.h"
#include "mapping_handler.h"
#include "multi_index_handler.h"
#include "src/spl/nurbs.h"

namespace iga {
template<int PARAMETRIC_DIMENSIONALITY>
class BasisFunctionHandler {
 public:
  explicit BasisFunctionHandler(std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>> spl) : spline_(std::move(spl)) {
    mapping_handler_ = std::make_shared<iga::MappingHandler<PARAMETRIC_DIMENSIONALITY>>(spline_);
    element_generator_ = std::make_shared<iga::elm::ElementGenerator<PARAMETRIC_DIMENSIONALITY>>(spline_);
  }

  std::vector<iga::elm::ElementIntegrationPoint<PARAMETRIC_DIMENSIONALITY>> EvaluateAllElementNonZeroNURBSBasisFunctions(
      int element_number, const iga::itg::IntegrationRule &rule) const {
    std::vector<iga::elm::ElementIntegrationPoint<PARAMETRIC_DIMENSIONALITY>> element_integration_points;
    std::array<std::vector<iga::itg::IntegrationPoint>, PARAMETRIC_DIMENSIONALITY> itg_pnts{};
    std::array<int, PARAMETRIC_DIMENSIONALITY> num_itg_pnts{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      itg_pnts[i] = rule.GetIntegrationPoints();
      num_itg_pnts[i] = rule.GetNumberOfIntegrationPoints();
    }
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> mih(num_itg_pnts);
    while (true) {
      std::array<double, PARAMETRIC_DIMENSIONALITY> itg_pnt_coords{};
      double itg_pnt_weight = 1;
      for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
        itg_pnt_coords[i] = itg_pnts[i][mih[i]].GetCoordinate();
        itg_pnt_weight *= itg_pnts[i][mih[i]].GetWeight();
      }
      std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coords = mapping_handler_->Reference2ParameterSpace(
          element_number, itg_pnt_coords);
      element_integration_points.emplace_back(iga::elm::ElementIntegrationPoint<PARAMETRIC_DIMENSIONALITY>(
          EvaluateAllNonZeroNURBSBasisFunctions(param_coords), itg_pnt_weight,
          mapping_handler_->GetJacobianDeterminant(param_coords)));
      if (mih.Get1DIndex() == mih.Get1DLength() - 1) break;
      ++mih;
    }
    return element_integration_points;
  }

  std::vector<iga::elm::ElementIntegrationPoint<PARAMETRIC_DIMENSIONALITY>> EvaluateAllElementNonZeroNURBSBasisFunctionDerivatives(
      int element_number, const iga::itg::IntegrationRule &rule) const {
    std::vector<iga::elm::ElementIntegrationPoint<PARAMETRIC_DIMENSIONALITY>> element_integration_points;
    std::array<std::vector<iga::itg::IntegrationPoint>, PARAMETRIC_DIMENSIONALITY> itg_pnts{};
    std::array<int, PARAMETRIC_DIMENSIONALITY> num_itg_pnts{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      itg_pnts[i] = rule.GetIntegrationPoints();
      num_itg_pnts[i] = rule.GetNumberOfIntegrationPoints();
    }
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> mih(num_itg_pnts);
    while (true) {
      std::array<double, PARAMETRIC_DIMENSIONALITY> itg_pnt_coords{};
      double itg_pnt_weight = 1;
      for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
        itg_pnt_coords[i] = itg_pnts[i][mih[i]].GetCoordinate();
        itg_pnt_weight *= itg_pnts[i][mih[i]].GetWeight();
      }
      std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coords = mapping_handler_->Reference2ParameterSpace(
          element_number, itg_pnt_coords);
      element_integration_points.emplace_back(iga::elm::ElementIntegrationPoint<PARAMETRIC_DIMENSIONALITY>(
          EvaluateAllNonZeroNURBSBasisFunctionDerivatives(param_coords), itg_pnt_weight,
          mapping_handler_->GetJacobianDeterminant(param_coords)));
      if (mih.Get1DIndex() == mih.Get1DLength() - 1) break;
      ++mih;
    }
    return element_integration_points;
  }

  std::vector<iga::elm::ElementIntegrationPoint<PARAMETRIC_DIMENSIONALITY>> EvaluateAllElementNonZeroNURBSBafDerivativesPhysical(
      int element_number, const iga::itg::IntegrationRule &rule) const {
    std::vector<iga::elm::ElementIntegrationPoint<PARAMETRIC_DIMENSIONALITY>> element_integration_points;
    std::array<std::vector<iga::itg::IntegrationPoint>, PARAMETRIC_DIMENSIONALITY> itg_pnts{};
    std::array<int, PARAMETRIC_DIMENSIONALITY> num_itg_pnts{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      itg_pnts[i] = rule.GetIntegrationPoints();
      num_itg_pnts[i] = rule.GetNumberOfIntegrationPoints();
    }
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> mih(num_itg_pnts);
    while (true) {
      std::array<double, PARAMETRIC_DIMENSIONALITY> itg_pnt_coords{};
      double itg_pnt_weight = 1;
      for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
        itg_pnt_coords[i] = itg_pnts[i][mih[i]].GetCoordinate();
        itg_pnt_weight *= itg_pnts[i][mih[i]].GetWeight();
      }
      std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coords = mapping_handler_->Reference2ParameterSpace(
          element_number, itg_pnt_coords);
      element_integration_points.emplace_back(iga::elm::ElementIntegrationPoint<PARAMETRIC_DIMENSIONALITY>(
          EvaluateAllNonZeroNURBSBafDerivativesPhysical(param_coords), itg_pnt_weight,
          mapping_handler_->GetJacobianDeterminant(param_coords)));
      if (mih.Get1DIndex() == mih.Get1DLength() - 1) break;
      ++mih;
    }
    return element_integration_points;
  }

 private:
  std::vector<double> EvaluateAllNonZeroNURBSBasisFunctions(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord) const {
    std::array<std::vector<double>, PARAMETRIC_DIMENSIONALITY> basis_functions{};
    std::array<int, PARAMETRIC_DIMENSIONALITY> num_baf{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      basis_functions[i] = spline_->EvaluateAllNonZeroBasisFunctions(i, param_coord[i]);
      num_baf[i] = basis_functions[i].size();
    }
    std::vector<double> nurbs_basis_functions;
    double sum = 0;
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> mih(num_baf);
    while (true) {
      double temp = 1;
      for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
        temp *= basis_functions[i][mih[i]];
      }
      nurbs_basis_functions.emplace_back(temp * GetWeight(param_coord, mih.Get1DIndex()));
      sum += temp * GetWeight(param_coord, mih.Get1DIndex());
      if (mih.Get1DIndex() == mih.Get1DLength() - 1) break;
      ++mih;
    }
    for (auto &nbaf : nurbs_basis_functions) {
      nbaf = nbaf / sum;
    }
    return nurbs_basis_functions;
  }

  std::array<std::vector<double>, PARAMETRIC_DIMENSIONALITY> EvaluateAllNonZeroNURBSBasisFunctionDerivatives(
      std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord) const {
    std::array<std::vector<double>, PARAMETRIC_DIMENSIONALITY> basis_functions{};
    std::array<std::vector<double>, PARAMETRIC_DIMENSIONALITY> basis_function_derivatives{};
    std::array<int, PARAMETRIC_DIMENSIONALITY> num_baf{};
    std::array<int, PARAMETRIC_DIMENSIONALITY> num_baf_ders{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      basis_functions[i] = spline_->EvaluateAllNonZeroBasisFunctions(i, param_coord[i]);
      basis_function_derivatives[i] = spline_->EvaluateAllNonZeroBasisFunctionDerivatives(i, param_coord[i], 1);
      num_baf[i] = basis_functions[i].size();
      num_baf_ders[i] = basis_function_derivatives[i].size();
    }
    std::vector<double> nurbs_basis_functions;
    std::array<std::vector<double>, PARAMETRIC_DIMENSIONALITY> nurbs_basis_function_derivatives;
    double sum_baf = 0;
    std::array<double, PARAMETRIC_DIMENSIONALITY> sum_baf_ders{};
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> mih(num_baf);
    while (true) {
      double temp = 1;
      std::array<double, PARAMETRIC_DIMENSIONALITY> temp_ders{};
      temp_ders.fill(1.0);
      for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
        temp *= basis_functions[i][mih[i]];
        for (int j = 0; j < PARAMETRIC_DIMENSIONALITY; ++j) {
          if (i == j) temp_ders[i] *= basis_function_derivatives[j][mih[j]];
          if (i != j) temp_ders[i] *= basis_functions[j][mih[j]];
        }
      }
      nurbs_basis_functions.emplace_back(temp * GetWeight(param_coord, mih.Get1DIndex()));
      sum_baf += temp * GetWeight(param_coord, mih.Get1DIndex());
      for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
        nurbs_basis_function_derivatives[i].emplace_back(temp_ders[i] * GetWeight(param_coord, mih.Get1DIndex()));
        sum_baf_ders[i] += temp_ders[i] * GetWeight(param_coord, mih.Get1DIndex());
      }
      if (mih.Get1DIndex() == mih.Get1DLength() - 1) break;
      ++mih;
    }
    for (uint64_t i = 0; i < basis_functions.size(); ++i) {
      nurbs_basis_functions[i] = nurbs_basis_functions[i] / sum_baf;
      for (int j = 0; j < PARAMETRIC_DIMENSIONALITY; ++j) {
        nurbs_basis_function_derivatives[j][i] = (nurbs_basis_function_derivatives[j][i] * sum_baf -
            nurbs_basis_functions[i] * sum_baf_ders[j]) / pow(sum_baf, 2);
      }
    }
    return nurbs_basis_function_derivatives;
  }

  std::array<std::vector<double>, PARAMETRIC_DIMENSIONALITY> EvaluateAllNonZeroNURBSBafDerivativesPhysical(
      std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord) const {
    std::array<std::vector<double>, PARAMETRIC_DIMENSIONALITY> dr_dx;
    std::array<std::vector<double>, PARAMETRIC_DIMENSIONALITY> dr_dxi = EvaluateAllNonZeroNURBSBasisFunctionDerivatives(param_coord);
    arma::dmat dxi_dx = mapping_handler_->GetDxiDx(param_coord);
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      for (uint64_t j = 0; j < dr_dxi[i].size(); ++j) {
        double temp = 0;
        for (int k = 0; k < PARAMETRIC_DIMENSIONALITY; ++k) {
          temp += dr_dxi[k][j] * dxi_dx(static_cast<uint64_t>(k), static_cast<uint64_t>(i));
        }
        dr_dx[i].emplace_back(temp);
      }
    }
    return dr_dx;
  }

  double GetWeight(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord, int local_index) const {
    iga::ConnectivityHandler<PARAMETRIC_DIMENSIONALITY> connectivity_handler(spline_);
    return spline_->GetWeight({connectivity_handler.GetGlobalIndex(element_generator_->GetElementNumberAtParametricCoordinate(
        param_coord), local_index) - 1});
  }

  std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>> spline_;
  std::shared_ptr<iga::MappingHandler<PARAMETRIC_DIMENSIONALITY>> mapping_handler_;
  std::shared_ptr<iga::elm::ElementGenerator<PARAMETRIC_DIMENSIONALITY>> element_generator_;
};
}  // namespace iga

#endif  // SRC_IGA_BASIS_FUNCTION_HANDLER_H_
