/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <armadillo>

#include "bdf_handler.h"
#include "four_point_gauss_legendre.h"
#include "gmock/gmock.h"
#include "linear_equation_assembler.h"
#include "nurbs.h"
#include "solution_vtk_writer.h"

using testing::Test;

class ALine : public Test {
 public:
  std::array<baf::KnotVector, 1> knot_vector =
      {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.4}, ParamCoord{0.5},
                        ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})};
  std::array<Degree, 1> degree = {Degree{3}};
  std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1};
  std::vector<baf::ControlPoint> control_points = {
      baf::ControlPoint(std::vector<double>({0.0})),
      baf::ControlPoint(std::vector<double>({1.0})),
      baf::ControlPoint(std::vector<double>({2.0})),
      baf::ControlPoint(std::vector<double>({3.0})),
      baf::ControlPoint(std::vector<double>({4.0})),
      baf::ControlPoint(std::vector<double>({5.0})),
      baf::ControlPoint(std::vector<double>({6.0}))};
  std::array<std::shared_ptr<baf::KnotVector>, 1> kv_ptr = {std::make_shared<baf::KnotVector>(knot_vector[0])};
  std::shared_ptr<spl::NURBS<1>> nurbs_ = std::make_shared<spl::NURBS<1>>(kv_ptr, degree, control_points, weights);

  iga::BDFHandler<1> bdf_handler = iga::BDFHandler<1>(nurbs_);
  iga::LinearEquationAssembler<1> linear_equation_assembler = iga::LinearEquationAssembler<1>(nurbs_);
  iga::ElementIntegralCalculator<1> elm_itg_calc = iga::ElementIntegralCalculator<1>(nurbs_);
  int n = nurbs_->GetNumberOfControlPoints();
  std::shared_ptr<arma::dmat> matA1 = std::make_shared<arma::dmat>(n, n, arma::fill::zeros);
  std::shared_ptr<arma::dmat> matA2 = std::make_shared<arma::dmat>(n, n, arma::fill::zeros);
  std::shared_ptr<arma::dvec> vecB1 = std::make_shared<arma::dvec>(n, arma::fill::zeros);
  std::shared_ptr<arma::dvec> vecB2 = std::make_shared<arma::dvec>(n, arma::fill::zeros);
  arma::dvec solution = arma::dvec(static_cast<uint64_t>(n), arma::fill::zeros);
  std::shared_ptr<arma::dvec> srcCp = std::make_shared<arma::dvec>(n, arma::fill::ones);
  iga::itg::IntegrationRule rule = iga::itg::FourPointGaussLegendre();
  iga::SolutionVTKWriter<1> solution_vtk_writer;
};

/*TEST_F(ALine, Test) { // NOLINT
  std::shared_ptr<arma::dvec> u_prev = std::make_shared<arma::dvec>(n, arma::fill::zeros);
  linear_equation_assembler.GetLeftSide(rule, matA1, elm_itg_calc);
  bdf_handler.GetMatrix(rule, matA2, 0.5);
  linear_equation_assembler.GetRightSide(rule, vecB1, elm_itg_calc, srcCp);
  (*matA1) = (*matA1) + (*matA2);
  solution_vtk_writer.WriteSolutionToVTK(nurbs_, solution, {{10}},
                                         "/Users/christophsusen/Desktop/line/solution_0.vtk");
  for (int i = 1; i < 1000; ++i) {
    (*vecB2) = (*vecB1) + (*matA2) * (*u_prev);
    linear_equation_assembler.SetZeroBC(matA1, vecB2);
    solution = arma::solve(*matA1, *vecB2);
    (*u_prev) = solution;
    solution_vtk_writer.WriteSolutionToVTK(nurbs_, solution, {{10}},
        "/Users/christophsusen/Desktop/line/solution_" + std::to_string(i) + ".vtk");
  }
}*/