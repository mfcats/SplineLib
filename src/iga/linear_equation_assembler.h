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
#include <vector>

#include "element_integral_calculator.h"
#include "element_generator.h"
#include "integration_rule.h"
#include "multi_index_handler.h"
#include "nurbs.h"

namespace iga {
template<int DIM>
class LinearEquationAssembler {
 public:
  explicit LinearEquationAssembler(std::shared_ptr<spl::NURBS<DIM>> spl) : spline_(std::move(spl)) {
    elm_gen_ = std::make_shared<iga::elm::ElementGenerator<DIM>>(spline_);
  }

  void GetLeftSide(const iga::itg::IntegrationRule &rule, const std::shared_ptr<arma::dmat> &matA,
      const iga::ElementIntegralCalculator<DIM> &elm_itg_calc) const {
    for (int e = 0; e < elm_gen_->GetNumberOfElements(); ++e) {
      elm_itg_calc.GetLaplaceElementIntegral(e, rule, matA);
    }
  }

  void GetRightSide(const iga::itg::IntegrationRule &rule, const std::shared_ptr<arma::dvec> &vecB,
      const iga::ElementIntegralCalculator<DIM> &elm_itg_calc, const std::shared_ptr<arma::dvec> &srcCp) const {
    for (int e = 0; e < elm_gen_->GetNumberOfElements(); ++e) {
      elm_itg_calc.GetLaplaceElementIntegral(e, rule, vecB, srcCp);
    }
  }

  void GetRightSideNeumann(const iga::itg::IntegrationRule &rule, const std::shared_ptr<arma::dvec> &vecB,
      const std::array<std::array<std::shared_ptr<arma::dvec>, 2>, DIM> &NeumannCp) const {
    std::array<std::array<std::shared_ptr<spl::NURBS<DIM - 1>>, 2>, DIM> boundary_splines = GetBoundarySplines();
    std::array<std::array<std::vector<int>, 2>, DIM> boundary_spline_connectivity = GetBoundarySplineConnectivity();
    for (int i = 0; i < boundary_splines.size(); ++i) {
      for (int j = 0; j < boundary_splines[i].size(); ++j) {
        iga::elm::ElementGenerator<DIM - 1> elm_gen(boundary_splines[i][j]);
        iga::BasisFunctionHandler<DIM - 1> baf_handler(boundary_splines[i][j]);
        iga::ConnectivityHandler<DIM - 1> connectivity_handler(boundary_splines[i][j]);
        for (int e = 0; e < elm_gen.GetNumberOfElements(); ++e) {
          std::vector<iga::elm::ElementIntegrationPoint<DIM - 1>> elm_intgr_pnts =
              baf_handler.EvaluateAllElementNonZeroNURBSBasisFunctions(e, rule);
          for (auto &p : elm_intgr_pnts) {
            double bc_int_pnt = 0;
            for (int l = 0; l < p.GetNumberOfNonZeroBasisFunctions(); ++l) {
              bc_int_pnt += p.GetBasisFunctionValue(l) *
                  (*(NeumannCp[i][j]))(static_cast<uint64_t>(connectivity_handler.GetGlobalIndex(e, l) - 1));
            }
            for (int k = 0; k < p.GetNumberOfNonZeroBasisFunctions(); ++k) {
              double temp = p.GetBasisFunctionValue(k) * bc_int_pnt * p.GetWeight() * p.GetJacobianDeterminant();
              (*vecB)(static_cast<uint64_t>(
                          boundary_spline_connectivity[i][j][connectivity_handler.GetGlobalIndex(e, k) - 1])) += temp;
            }
          }
        }
      }
    }
  }

  std::array<std::array<std::shared_ptr<spl::NURBS<DIM - 1>>, 2>, DIM> GetBoundarySplines() const {
    std::array<int, DIM> points_per_dir = spline_->GetPointsPerDirection();
    std::array<std::array<std::shared_ptr<spl::NURBS<DIM - 1>>, 2>, DIM> boundary_splines;
    std::array<std::array<std::vector<int>, 2>, DIM> boundary_spline_connectivity = GetBoundarySplineConnectivity();
    for (int i = 0; i < DIM; ++i) {
      std::array<Degree, DIM - 1> degree;
      std::array<std::shared_ptr<baf::KnotVector>, DIM - 1> kv_ptr;
      int m = 0;
      for (int j = 0; j < DIM; ++j) {
        if (j != i) {
          kv_ptr[m] = spline_->GetKnotVector(m);
          degree[m] = spline_->GetDegree(m);
          ++m;
        }
      }
      std::array<std::vector<baf::ControlPoint>, 2> control_points;
      std::array<std::vector<double>, 2> weights;
      util::MultiIndexHandler<DIM> mih(points_per_dir);
      for (int k = 0; k < boundary_spline_connectivity[i].size(); ++k) {
        for (int n = 0; n < boundary_spline_connectivity[i][k].size(); ++n) {
          mih.Set1DIndex(boundary_spline_connectivity[i][k][n]);
          control_points[k].emplace_back(spline_->GetControlPoint(mih.GetIndices()));
          weights[k].emplace_back(spline_->GetWeight(mih.GetIndices()));
        }
        boundary_splines[i][k] = std::make_shared<spl::NURBS<DIM - 1>>(kv_ptr, degree, control_points[k], weights[k]);
      }
    }
    return boundary_splines;
  }

  std::array<std::array<std::vector<int>, 2>, DIM> GetBoundarySplineConnectivity() const {
    std::array<int, DIM> points_per_dir = spline_->GetPointsPerDirection();
    std::array<std::array<std::vector<int>, 2>, DIM> boundary_spl_connectivity;
    for (int i = 0; i < DIM; ++i) {
      util::MultiIndexHandler<DIM> mih(points_per_dir);
      for (int k = 0; k < mih.Get1DLength(); ++k) {
        if (mih[i] == 0) {
          boundary_spl_connectivity[i][0].emplace_back(mih.Get1DIndex());
        }
        if (mih[i] == points_per_dir[i] - 1) {
          boundary_spl_connectivity[i][1].emplace_back(mih.Get1DIndex());
        }
        ++mih;
      }
    }
    return boundary_spl_connectivity;
  }

  void SetZeroBC(const std::shared_ptr<arma::dmat> &matA, const std::shared_ptr<arma::dvec> &vecB) {
    util::MultiIndexHandler<DIM> mih(spline_->GetPointsPerDirection());
    while (true) {
      bool on_boundary = false;
      for (int i = 0; i < DIM; ++i) {
        if (!((mih[i] > 0) && (mih.GetDifferenceIndices()[i] > 0))) {
          on_boundary = true;
        }
      }
      if (on_boundary) {
        (*vecB)(static_cast<uint64_t>(mih.Get1DIndex())) = 0;
        (*matA).row(static_cast<uint64_t>(mih.Get1DIndex())).fill(0);
        (*matA)(static_cast<uint64_t>(mih.Get1DIndex()), static_cast<uint64_t>(mih.Get1DIndex())) = 1;
      }
      if (mih.Get1DIndex() == mih.Get1DLength() - 1) break;
      ++mih;
    }
  }

  // Only used for test case in test/solution_vtk_writer_examples.cc which is currently commented out.
  /*
  void SetLinearBC(const std::shared_ptr<arma::dmat> &matA, const std::shared_ptr<arma::dvec> &vecB) {
    uint64_t l = 0;
    uint64_t k = 0;
    int n = spline_->GetPointsPerDirection()[0];
    int m = spline_->GetPointsPerDirection()[1];
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < n; ++i) {
        if (i == 0) {
          (*vecB)(l) = 0;
          (*matA).row(l).fill(0);
          (*matA)(l, l) = 1;
        } else if (i == n - 1) {
          (*vecB)(l) = 1;
          (*matA).row(l).fill(0);
          (*matA)(l, l) = 1;
        } else if ((j == 0) || (j == m - 1)) {
          (*vecB)(l) = 0 + ((1.0 - 0.0) / (2.0 - 0.0)) * (spline_->GetControlPoints()[k] - 0.0);
          (*matA).row(l).fill(0);
          (*matA)(l, l) = 1;
        }
        l += 1;
        k += 3;
      }
    }
  }
  */

 private:
  std::shared_ptr<spl::NURBS<DIM>> spline_;
  std::shared_ptr<iga::elm::ElementGenerator<DIM>> elm_gen_;
};
}  // namespace iga

#endif  // SRC_IGA_LINEAR_EQUATION_ASSEMBLER_H_
