/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include <armadillo>

#include "gmock/gmock.h"
#include "solution_vtk_writer.h"
#include "test_spline.h"

// only needed for tests that are currently commented out
// #include "four_point_gauss_legendre.h"

TEST_F(AnIGATestSpline, TestSolutionVTKWriter) {  // NOLINT
  linear_equation_assembler.GetLeftSide(rule, matA, elm_itg_calc);
  linear_equation_assembler.GetRightSide(rule, vecB, elm_itg_calc, srcCp);
  linear_equation_assembler.SetZeroBC(matA, vecB);
  arma::dvec solution = arma::solve(*matA, *vecB);
  iga::SolutionVTKWriter<2> solution_vtk_writer;
  solution_vtk_writer.WriteSolutionToVTK(nurbs_, solution, {{10, 10}}, "solution.vtk");
  remove("solution.vtk");
}

// The tests below compute the solution of Laplace's equation on a line and inside a cube.

/*class ALine : public Test {
 public:
  std::array<baf::KnotVector, 1> knot_vector =
      {baf::KnotVector({ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.4}, ParametricCoordinate{0.5},
                        ParametricCoordinate{0.6}, ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{1}})};
  std::array<Degree, 1> degree = {Degree{3}};
  std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1};
  std::vector<spl::ControlPoint> control_points = {
      spl::ControlPoint(std::vector<double>({0.0})),
      spl::ControlPoint(std::vector<double>({1.0})),
      spl::ControlPoint(std::vector<double>({2.0})),
      spl::ControlPoint(std::vector<double>({3.0})),
      spl::ControlPoint(std::vector<double>({4.0})),
      spl::ControlPoint(std::vector<double>({5.0})),
      spl::ControlPoint(std::vector<double>({6.0}))};
  std::array<std::shared_ptr<baf::KnotVector>, 1> kv_ptr = {std::make_shared<baf::KnotVector>(knot_vector[0])};
  std::shared_ptr<spl::NURBS<1>> nurbs_ = std::make_shared<spl::NURBS<1>>(kv_ptr, degree, control_points, weights);

  iga::LinearEquationAssembler<1> linear_equation_assembler = iga::LinearEquationAssembler<1>(nurbs_);
  iga::ElementIntegralCalculator<1> elm_itg_calc = iga::ElementIntegralCalculator<1>(nurbs_);
  int n = nurbs_->GetNumberOfControlPoints();
  std::shared_ptr<arma::dmat> matA = std::make_shared<arma::dmat>(n, n, arma::fill::zeros);
  std::shared_ptr<arma::dvec> vecB = std::make_shared<arma::dvec>(n, arma::fill::zeros);
  std::shared_ptr<arma::dvec> srcCp = std::make_shared<arma::dvec>(n, arma::fill::ones);
  iga::itg::IntegrationRule rule = iga::itg::FourPointGaussLegendre();
};

TEST_F(ALine, Test) {  // NOLINT
  linear_equation_assembler.GetLeftSide(rule, matA, elm_itg_calc);
  linear_equation_assembler.GetRightSide(rule, vecB, elm_itg_calc, srcCp);
  linear_equation_assembler.SetZeroBC(matA, vecB);
  arma::dvec solution = arma::solve(*matA, *vecB);
  iga::SolutionVTKWriter<1> solution_vtk_writer;
  solution_vtk_writer.WriteSolutionToVTK(nurbs_, solution, {{10}}, "/Users/christophsusen/Desktop/solution.vtk");
}

class ACube : public Test {
 public:
  std::array<baf::KnotVector, 3> knot_vector =
      {baf::KnotVector({ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.5}, ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{1}}),
       baf::KnotVector({ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.5}, ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{1}}),
       baf::KnotVector({ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0}, ParametricCoordinate{0.5}, ParametricCoordinate{1}, ParametricCoordinate{1}, ParametricCoordinate{1}})};

  std::array<Degree, 3> degree = {Degree{2}, Degree{2}, Degree{2}};

  std::vector<double> weights = {1,1,1,1,1,1,1,1,
                                 1,1,1,1,1,1,1,1,
                                 1,1,1,1,1,1,1,1,
                                 1,1,1,1,1,1,1,1,
                                 1,1,1,1,1,1,1,1,
                                 1,1,1,1,1,1,1,1,
                                 1,1,1,1,1,1,1,1,
                                 1,1,1,1,1,1,1,1};

  std::vector<spl::ControlPoint> control_points = {
      spl::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
      spl::ControlPoint(std::vector<double>({0.33, 0.0, 0.0})),
      spl::ControlPoint(std::vector<double>({0.66, 0.0, 0.0})),
      spl::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
      spl::ControlPoint(std::vector<double>({0.0, 0.33, 0.0})),
      spl::ControlPoint(std::vector<double>({0.33, 0.33, 0.0})),
      spl::ControlPoint(std::vector<double>({0.66, 0.33, 0.0})),
      spl::ControlPoint(std::vector<double>({1.0, 0.33, 0.0})),
      spl::ControlPoint(std::vector<double>({0.0, 0.66, 0.0})),
      spl::ControlPoint(std::vector<double>({0.33, 0.66, 0.0})),
      spl::ControlPoint(std::vector<double>({0.66, 0.66, 0.0})),
      spl::ControlPoint(std::vector<double>({1.0, 0.66, 0.0})),
      spl::ControlPoint(std::vector<double>({0.0, 1.0, 0.0})),
      spl::ControlPoint(std::vector<double>({0.33, 1.0, 0.0})),
      spl::ControlPoint(std::vector<double>({0.66, 1.0, 0.0})),
      spl::ControlPoint(std::vector<double>({1.0, 1.0, 0.0})),

      spl::ControlPoint(std::vector<double>({0.0, 0.0, 0.33})),
      spl::ControlPoint(std::vector<double>({0.33, 0.0, 0.33})),
      spl::ControlPoint(std::vector<double>({0.66, 0.0, 0.33})),
      spl::ControlPoint(std::vector<double>({1.0, 0.0, 0.33})),
      spl::ControlPoint(std::vector<double>({0.0, 0.33, 0.33})),
      spl::ControlPoint(std::vector<double>({0.33, 0.33, 0.33})),
      spl::ControlPoint(std::vector<double>({0.66, 0.33, 0.33})),
      spl::ControlPoint(std::vector<double>({1.0, 0.33, 0.33})),
      spl::ControlPoint(std::vector<double>({0.0, 0.66, 0.33})),
      spl::ControlPoint(std::vector<double>({0.33, 0.66, 0.33})),
      spl::ControlPoint(std::vector<double>({0.66, 0.66, 0.33})),
      spl::ControlPoint(std::vector<double>({1.0, 0.66, 0.33})),
      spl::ControlPoint(std::vector<double>({0.0, 1.0, 0.33})),
      spl::ControlPoint(std::vector<double>({0.33, 1.0, 0.33})),
      spl::ControlPoint(std::vector<double>({0.66, 1.0, 0.33})),
      spl::ControlPoint(std::vector<double>({1.0, 1.0, 0.33})),

      spl::ControlPoint(std::vector<double>({0.0, 0.0, 0.66})),
      spl::ControlPoint(std::vector<double>({0.33, 0.0, 0.66})),
      spl::ControlPoint(std::vector<double>({0.66, 0.0, 0.66})),
      spl::ControlPoint(std::vector<double>({1.0, 0.0, 0.66})),
      spl::ControlPoint(std::vector<double>({0.0, 0.33, 0.66})),
      spl::ControlPoint(std::vector<double>({0.33, 0.33, 0.66})),
      spl::ControlPoint(std::vector<double>({0.66, 0.33, 0.66})),
      spl::ControlPoint(std::vector<double>({1.0, 0.33, 0.66})),
      spl::ControlPoint(std::vector<double>({0.0, 0.66, 0.66})),
      spl::ControlPoint(std::vector<double>({0.33, 0.66, 0.66})),
      spl::ControlPoint(std::vector<double>({0.66, 0.66, 0.66})),
      spl::ControlPoint(std::vector<double>({1.0, 0.66, 0.66})),
      spl::ControlPoint(std::vector<double>({0.0, 1.0, 0.66})),
      spl::ControlPoint(std::vector<double>({0.33, 1.0, 0.66})),
      spl::ControlPoint(std::vector<double>({0.66, 1.0, 0.66})),
      spl::ControlPoint(std::vector<double>({1.0, 1.0, 0.66})),

      spl::ControlPoint(std::vector<double>({0.0, 0.0, 1.0})),
      spl::ControlPoint(std::vector<double>({0.33, 0.0, 1.0})),
      spl::ControlPoint(std::vector<double>({0.66, 0.0, 1.0})),
      spl::ControlPoint(std::vector<double>({1.0, 0.0, 1.0})),
      spl::ControlPoint(std::vector<double>({0.0, 0.33, 1.0})),
      spl::ControlPoint(std::vector<double>({0.33, 0.33, 1.0})),
      spl::ControlPoint(std::vector<double>({0.66, 0.33, 1.0})),
      spl::ControlPoint(std::vector<double>({1.0, 0.33, 1.0})),
      spl::ControlPoint(std::vector<double>({0.0, 0.66, 1.0})),
      spl::ControlPoint(std::vector<double>({0.33, 0.66, 1.0})),
      spl::ControlPoint(std::vector<double>({0.66, 0.66, 1.0})),
      spl::ControlPoint(std::vector<double>({1.0, 0.66, 1.0})),
      spl::ControlPoint(std::vector<double>({0.0, 1.0, 1.0})),
      spl::ControlPoint(std::vector<double>({0.33, 1.0, 1.0})),
      spl::ControlPoint(std::vector<double>({0.66, 1.0, 1.0})),
      spl::ControlPoint(std::vector<double>({1.0, 1.0, 1.0}))
  };

  std::array<std::shared_ptr<baf::KnotVector>, 3> kv_ptr = {std::make_shared<baf::KnotVector>(knot_vector[0]),
                                                            std::make_shared<baf::KnotVector>(knot_vector[1]),
                                                            std::make_shared<baf::KnotVector>(knot_vector[2])};
  std::shared_ptr<spl::NURBS<3>> nurbs_ = std::make_shared<spl::NURBS<3>>(kv_ptr, degree, control_points, weights);

  iga::LinearEquationAssembler<3> linear_equation_assembler = iga::LinearEquationAssembler<3>(nurbs_);
  iga::ElementIntegralCalculator<3> elm_itg_calc = iga::ElementIntegralCalculator<3>(nurbs_);
  int n = nurbs_->GetNumberOfControlPoints();
  std::shared_ptr<arma::dmat> matA = std::make_shared<arma::dmat>(n, n, arma::fill::zeros);
  std::shared_ptr<arma::dvec> vecB = std::make_shared<arma::dvec>(n, arma::fill::zeros);
  std::shared_ptr<arma::dvec> srcCp = std::make_shared<arma::dvec>(n, arma::fill::ones);
  iga::itg::IntegrationRule rule = iga::itg::TwoPointGaussLegendre();
};

TEST_F(ACube, Test) {  // NOLINT
  linear_equation_assembler.GetLeftSide(rule, matA, elm_itg_calc);
  linear_equation_assembler.GetRightSide(rule, vecB, elm_itg_calc, srcCp);
  linear_equation_assembler.SetZeroBC(matA, vecB);
  arma::dvec solution = arma::solve(*matA, *vecB);
  iga::SolutionVTKWriter<3> solution_vtk_writer;
  solution_vtk_writer.WriteSolutionToVTK(nurbs_, solution, {{10,10,10}}, "/Users/christophsusen/Desktop/solution.vtk");
}*/
