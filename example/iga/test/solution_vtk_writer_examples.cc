/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

/*
 * The following two tests are commented out to reduce the time that is needed to run all travis ci tests.
 * All code that is used by these two tests is already tested with less complex test cases.
 * The first test shows how to refine the l-shaped nurbs before the solution is computed.
 * The second test computes the solution of the laplace equation on a rectangle with boundary condition one on one side
 * and boundary condition zero on the other side.
 * The results are written to a vtk file and can be viewed using paraview.
*/

/*
#include <armadillo>

#include "gmock/gmock.h"
#include "solution_vtk_writer.h"
#include "test_spline.h"
#include "four_point_gauss_legendre.h"

TEST_F(AnIGATestSpline, TestSolutionVTKWriterRefined) { // NOLINT
  std::shared_ptr<spl::NURBS<2>> nurbs_refined = nurbs_;
  std::vector<double> knots_to_add = {0.48, 0.49, 0.51, 0.52};
  for (auto &knot : knots_to_add) {
    nurbs_refined->InsertKnot(ParametricCoordinate(knot), 0);
    nurbs_refined->InsertKnot(ParametricCoordinate(knot), 1);
  }
  iga::LinearEquationAssembler linear_equation_assembler_ref(nurbs_refined);
  int n_ref = nurbs_refined->GetNumberOfControlPoints();
  std::shared_ptr<arma::dmat> matA_ref = std::make_shared<arma::dmat>(n_ref, n_ref, arma::fill::zeros);
  std::shared_ptr<arma::dvec> vecB_ref = std::make_shared<arma::dvec>(n_ref, arma::fill::zeros);
  std::shared_ptr<arma::dvec> srcCp_ref = std::make_shared<arma::dvec>(n_ref, arma::fill::ones);
  iga::itg::IntegrationRule rule_ref = iga::itg::FourPointGaussLegendre();
  iga::ElementIntegralCalculator elm_itg_calc_ref = iga::ElementIntegralCalculator(nurbs_refined);
  linear_equation_assembler_ref.GetLeftSide(rule_ref, matA_ref, elm_itg_calc_ref);
  linear_equation_assembler_ref.GetRightSide(rule_ref, vecB_ref, elm_itg_calc_ref, srcCp_ref);
  linear_equation_assembler_ref.SetZeroBC(matA_ref, vecB_ref);
  arma::dvec solution_ref = arma::solve(*matA_ref, *vecB_ref);
  iga::SolutionVTKWriter solution_vtk_writer;
  solution_vtk_writer.WriteSolutionToVTK(nurbs_refined, solution_ref, {{30, 30}}, "solution_refined.vtk");
}

TEST_F(AnIGATestSpline2, TestSolutionVTKWriter) { // NOLINT
  linear_equation_assembler.GetLeftSide(rule, matA, elm_itg_calc);
  linear_equation_assembler.GetRightSide(rule, vecB, elm_itg_calc, srcCp);
  linear_equation_assembler.SetLinearBC(matA, vecB);
  arma::dvec solution = arma::solve(*matA, *vecB);
  iga::SolutionVTKWriter solution_vtk_writer;
  solution_vtk_writer.WriteSolutionToVTK(nurbs_, solution, {{30, 30}}, "solution_2.vtk");
}
*/
