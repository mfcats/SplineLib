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
#include "element_integral_calculator.h"
#include "four_point_gauss_legendre.h"
#include "gmock/gmock.h"
#include "linear_equation_assembler.h"
#include "nurbs.h"
#include "poisson_problem.h"
#include "solution_vtk_writer.h"

using testing::Test;
using testing::DoubleNear;

class ALine : public Test {
 public:
  std::array<baf::KnotVector, 1> knot_vector =
      {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.1}, ParamCoord{0.2},
                        ParamCoord{0.3}, ParamCoord{0.4}, ParamCoord{0.5}, ParamCoord{0.6}, ParamCoord{0.7},
                        ParamCoord{0.8}, ParamCoord{0.9}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})};
  std::array<Degree, 1> degree = {Degree{3}};
  std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<baf::ControlPoint> control_points = {
      baf::ControlPoint(std::vector<double>({0.0})),
      baf::ControlPoint(std::vector<double>({0.5})),
      baf::ControlPoint(std::vector<double>({1.0})),
      baf::ControlPoint(std::vector<double>({1.5})),
      baf::ControlPoint(std::vector<double>({2.0})),
      baf::ControlPoint(std::vector<double>({3.0})),
      baf::ControlPoint(std::vector<double>({4.0})),
      baf::ControlPoint(std::vector<double>({5.0})),
      baf::ControlPoint(std::vector<double>({6.0})),
      baf::ControlPoint(std::vector<double>({6.5})),
      baf::ControlPoint(std::vector<double>({7.0})),
      baf::ControlPoint(std::vector<double>({7.5})),
      baf::ControlPoint(std::vector<double>({8.0}))};
  std::array<std::shared_ptr<baf::KnotVector>, 1> kv_ptr = {std::make_shared<baf::KnotVector>(knot_vector[0])};
  std::shared_ptr<spl::NURBS<1>> nurbs_ = std::make_shared<spl::NURBS<1>>(kv_ptr, degree, control_points, weights);
};

TEST_F(ALine, Test) { // NOLINT
  iga::itg::IntegrationRule rule = iga::itg::FourPointGaussLegendre();

  iga::PoissonProblem<1> poisson_problem(nurbs_, rule);
  auto solutions = poisson_problem.GetUnsteadyStateSolution(0.5, 100);
  auto steady_solution = poisson_problem.GetSteadyStateSolution();
  for (uint64_t i = 0; i < steady_solution.size(); ++i) {
    ASSERT_THAT((*solutions[solutions.size() - 1])(i), DoubleNear(steady_solution(i), 0.0005));
  }

  // Code to write all time steps to vtk files.
  /*iga::SolutionVTKWriter<1> solution_vtk_writer;
  for (uint64_t i = 0; i < solutions.size(); ++i) {
  solution_vtk_writer.WriteSolutionToVTK(nurbs_, *(solutions[i]), {{30}},
      "/Users/christophsusen/Desktop/line/solution_" + std::to_string(i) + ".vtk");
}*/
}

/*class ASquarePlate : public Test {
 public:
  std::array<baf::KnotVector, 2> knot_vector =
      {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.4}, ParamCoord{0.5},
                        ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}),
       baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.4}, ParamCoord{0.5},
                        ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})};
  std::array<Degree, 2> degree = {Degree{3}, Degree{3}};
  std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1};
  std::vector<baf::ControlPoint> control_points = {
      baf::ControlPoint(std::vector<double>({0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({2.0, 0.0})),
      baf::ControlPoint(std::vector<double>({3.0, 0.0})),
      baf::ControlPoint(std::vector<double>({4.0, 0.0})),
      baf::ControlPoint(std::vector<double>({5.0, 0.0})),
      baf::ControlPoint(std::vector<double>({6.0, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, 1.0})),
      baf::ControlPoint(std::vector<double>({1.0, 1.0})),
      baf::ControlPoint(std::vector<double>({2.0, 1.0})),
      baf::ControlPoint(std::vector<double>({3.0, 1.0})),
      baf::ControlPoint(std::vector<double>({4.0, 1.0})),
      baf::ControlPoint(std::vector<double>({5.0, 1.0})),
      baf::ControlPoint(std::vector<double>({6.0, 1.0})),

      baf::ControlPoint(std::vector<double>({0.0, 2.0})),
      baf::ControlPoint(std::vector<double>({1.0, 2.0})),
      baf::ControlPoint(std::vector<double>({2.0, 2.0})),
      baf::ControlPoint(std::vector<double>({3.0, 2.0})),
      baf::ControlPoint(std::vector<double>({4.0, 2.0})),
      baf::ControlPoint(std::vector<double>({5.0, 2.0})),
      baf::ControlPoint(std::vector<double>({6.0, 2.0})),

      baf::ControlPoint(std::vector<double>({0.0, 3.0})),
      baf::ControlPoint(std::vector<double>({1.0, 3.0})),
      baf::ControlPoint(std::vector<double>({2.0, 3.0})),
      baf::ControlPoint(std::vector<double>({3.0, 3.0})),
      baf::ControlPoint(std::vector<double>({4.0, 3.0})),
      baf::ControlPoint(std::vector<double>({5.0, 3.0})),
      baf::ControlPoint(std::vector<double>({6.0, 3.0})),

      baf::ControlPoint(std::vector<double>({0.0, 4.0})),
      baf::ControlPoint(std::vector<double>({1.0, 4.0})),
      baf::ControlPoint(std::vector<double>({2.0, 4.0})),
      baf::ControlPoint(std::vector<double>({3.0, 4.0})),
      baf::ControlPoint(std::vector<double>({4.0, 4.0})),
      baf::ControlPoint(std::vector<double>({5.0, 4.0})),
      baf::ControlPoint(std::vector<double>({6.0, 4.0})),

      baf::ControlPoint(std::vector<double>({0.0, 5.0})),
      baf::ControlPoint(std::vector<double>({1.0, 5.0})),
      baf::ControlPoint(std::vector<double>({2.0, 5.0})),
      baf::ControlPoint(std::vector<double>({3.0, 5.0})),
      baf::ControlPoint(std::vector<double>({4.0, 5.0})),
      baf::ControlPoint(std::vector<double>({5.0, 5.0})),
      baf::ControlPoint(std::vector<double>({6.0, 5.0})),

      baf::ControlPoint(std::vector<double>({0.0, 6.0})),
      baf::ControlPoint(std::vector<double>({1.0, 6.0})),
      baf::ControlPoint(std::vector<double>({2.0, 6.0})),
      baf::ControlPoint(std::vector<double>({3.0, 6.0})),
      baf::ControlPoint(std::vector<double>({4.0, 6.0})),
      baf::ControlPoint(std::vector<double>({5.0, 6.0})),
      baf::ControlPoint(std::vector<double>({6.0, 6.0}))};
  std::array<std::shared_ptr<baf::KnotVector>, 2> kv_ptr = {std::make_shared<baf::KnotVector>(knot_vector[0]),
                                                            std::make_shared<baf::KnotVector>(knot_vector[1])};
  std::shared_ptr<spl::NURBS<2>> nurbs_ = std::make_shared<spl::NURBS<2>>(kv_ptr, degree, control_points, weights);
};

TEST_F(ASquarePlate, Test2) { // NOLINT
  iga::itg::IntegrationRule rule = iga::itg::FourPointGaussLegendre();
  iga::SolutionVTKWriter<2> solution_vtk_writer;
  iga::PoissonProblem<2> poisson_problem(nurbs_, rule);
  auto solutions = poisson_problem.GetUnsteadyStateSolution(0.5, 25);
  for (uint64_t i = 0; i < solutions.size(); ++i) {
    solution_vtk_writer.WriteSolutionToVTK(nurbs_, *(solutions[i]), {{30, 30}},
        "/Users/christophsusen/Desktop/square/solution_" + std::to_string(i) + ".vtk");
  }
}*/

/*TEST_F(ASquarePlate, Test3) { // NOLINT
  io::VTKWriter vtk_writer;
  iga::LinearEquationAssembler<2> lin_eq_assem(nurbs_);
  std::array<std::shared_ptr<spl::NURBS<1>>, 4> boundary_splines = lin_eq_assem.GetBoundarySplines();
  for (uint64_t i = 0; i < boundary_splines.size(); ++i) {
    auto spl = std::make_any<std::shared_ptr<spl::NURBS<1>>>(boundary_splines[i]);
    vtk_writer.WriteFile({spl}, "/Users/christophsusen/Desktop/boundary_spline/spline_" + std::to_string(i) + ".vtk",
                         {{30}, {30, 30}, {30, 30, 30}});
  }
}*/
