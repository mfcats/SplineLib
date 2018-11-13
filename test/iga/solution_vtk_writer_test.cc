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

#include "gmock/gmock.h"
#include "solution_vtk_writer.h"
#include "test_spline.h"

TEST_F(AnIGATestSpline, TestSolutionVTKWriter) { // NOLINT
  linear_equation_assembler.GetLeftSide(rule, matA, elm_itg_calc);
  linear_equation_assembler.GetRightSide(rule, vecB, elm_itg_calc, srcCp);
  linear_equation_assembler.SetZeroBC(matA, vecB);
  arma::dvec solution = arma::solve(*matA, *vecB);
  iga::SolutionVTKWriter solution_vtk_writer;
  solution_vtk_writer.WriteSolutionToVTK(nurbs_, solution, {{10, 10}}, "solution.vtk");
  remove("solution.vtk");
}

TEST_F(AnIGATestSpline2, TestSolutionVTKWriter) { // NOLINT
  linear_equation_assembler.GetLeftSide(rule, matA, elm_itg_calc);
  linear_equation_assembler.GetRightSide(rule, vecB, elm_itg_calc, srcCp);
  linear_equation_assembler.SetLinearBC(matA, vecB);
  arma::dvec solution = arma::solve(*matA, *vecB);
  iga::SolutionVTKWriter solution_vtk_writer;
  solution_vtk_writer.WriteSolutionToVTK(nurbs_, solution, {{30, 30}}, "solution_2.vtk");
  remove("solution_2.vtk");
}
