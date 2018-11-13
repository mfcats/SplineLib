/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef TEST_IGA_TEST_SPLINE_H_
#define TEST_IGA_TEST_SPLINE_H_

#include <vector>

#include "gmock/gmock.h"
#include "basis_function_handler.h"
#include "linear_equation_assembler.h"
#include "nurbs.h"
#include "two_point_gauss_legendre.h"
#include "four_point_gauss_legendre.h"

using testing::Test;

class AnIGATestSpline : public Test {
 public:
  std::array<baf::KnotVector, 2> knot_vector =
      {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.4}, ParamCoord{0.5},
                        ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}),
       baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.5}, ParamCoord{0.5},
                        ParamCoord{0.5}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})};

  std::array<Degree, 2> degree = {Degree{3}, Degree{3}};

  std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1};

  std::vector<baf::ControlPoint> control_points = {
      baf::ControlPoint(std::vector<double>({0.0, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.16, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.32, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.48, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.64, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.8, -1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -1.0, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.16, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.32, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.48, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.64, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({0.8, -0.66, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -0.66, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.16, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.32, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.48, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.64, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({0.8, -0.33, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, -0.33, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.16, 0.16, 0.0})),
      baf::ControlPoint(std::vector<double>({0.32, 0.32, 0.0})),
      baf::ControlPoint(std::vector<double>({0.48, 0.48, 0.0})),
      baf::ControlPoint(std::vector<double>({0.64, 0.64, 0.0})),
      baf::ControlPoint(std::vector<double>({0.8, 0.8, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0})),

      baf::ControlPoint(std::vector<double>({-0.33, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.16, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.32, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.48, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.64, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 0.8, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.33, 1.0, 0.0})),

      baf::ControlPoint(std::vector<double>({-0.66, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.16, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.32, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.48, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.64, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 0.8, 0.0})),
      baf::ControlPoint(std::vector<double>({-0.66, 1.0, 0.0})),

      baf::ControlPoint(std::vector<double>({-1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.16, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.32, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.48, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.64, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 0.8, 0.0})),
      baf::ControlPoint(std::vector<double>({-1.0, 1.0, 0.0}))
  };

  std::array<std::shared_ptr<baf::KnotVector>, 2> kv_ptr = {std::make_shared<baf::KnotVector>(knot_vector[0]),
                                                            std::make_shared<baf::KnotVector>(knot_vector[1])};
  std::shared_ptr<spl::NURBS<2>> nurbs_ = std::make_shared<spl::NURBS<2>>(kv_ptr, degree, control_points, weights);

  iga::LinearEquationAssembler linear_equation_assembler = iga::LinearEquationAssembler(nurbs_);
  iga::ElementIntegralCalculator elm_itg_calc = iga::ElementIntegralCalculator(nurbs_);
  int n = nurbs_->GetNumberOfControlPoints();
  std::shared_ptr<arma::dmat> matA = std::make_shared<arma::dmat>(n, n, arma::fill::zeros);
  std::shared_ptr<arma::dvec> vecB = std::make_shared<arma::dvec>(n, arma::fill::zeros);
  std::shared_ptr<arma::dvec> srcCp = std::make_shared<arma::dvec>(n, arma::fill::ones);
  iga::itg::IntegrationRule rule = iga::itg::TwoPointGaussLegendre();
  iga::BasisFunctionHandler basis_function_handler = iga::BasisFunctionHandler(nurbs_);
};

class AnIGATestSpline2 : public Test {
 public:
  std::array<baf::KnotVector, 2> knot_vector =
      {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.15}, ParamCoord{0.3},
                        ParamCoord{0.45}, ParamCoord{0.55}, ParamCoord{0.7}, ParamCoord{0.85}, ParamCoord{1},
                        ParamCoord{1}, ParamCoord{1}}),
       baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.15}, ParamCoord{0.3},
                        ParamCoord{0.45}, ParamCoord{0.55}, ParamCoord{0.7}, ParamCoord{0.85}, ParamCoord{1},
                        ParamCoord{1}, ParamCoord{1}})};

  std::array<Degree, 2> degree = {Degree{2}, Degree{2}};

  std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1};

  std::vector<baf::ControlPoint> control_points = {
      baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.25, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.5, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.75, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.25, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.5, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.75, 0.0, 0.0})),
      baf::ControlPoint(std::vector<double>({2.0, 0.0, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, 0.125, 0.0})),
      baf::ControlPoint(std::vector<double>({0.25, 0.125, 0.0})),
      baf::ControlPoint(std::vector<double>({0.5, 0.125, 0.0})),
      baf::ControlPoint(std::vector<double>({0.75, 0.125, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.125, 0.0})),
      baf::ControlPoint(std::vector<double>({1.25, 0.125, 0.0})),
      baf::ControlPoint(std::vector<double>({1.5, 0.125, 0.0})),
      baf::ControlPoint(std::vector<double>({1.75, 0.125, 0.0})),
      baf::ControlPoint(std::vector<double>({2.0, 0.125, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, 0.25, 0.0})),
      baf::ControlPoint(std::vector<double>({0.25, 0.25, 0.0})),
      baf::ControlPoint(std::vector<double>({0.5, 0.25, 0.0})),
      baf::ControlPoint(std::vector<double>({0.75, 0.25, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.25, 0.0})),
      baf::ControlPoint(std::vector<double>({1.25, 0.25, 0.0})),
      baf::ControlPoint(std::vector<double>({1.5, 0.25, 0.0})),
      baf::ControlPoint(std::vector<double>({1.75, 0.25, 0.0})),
      baf::ControlPoint(std::vector<double>({2.0, 0.25, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, 0.375, 0.0})),
      baf::ControlPoint(std::vector<double>({0.25, 0.375, 0.0})),
      baf::ControlPoint(std::vector<double>({0.5, 0.375, 0.0})),
      baf::ControlPoint(std::vector<double>({0.75, 0.375, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.375, 0.0})),
      baf::ControlPoint(std::vector<double>({1.25, 0.375, 0.0})),
      baf::ControlPoint(std::vector<double>({1.5, 0.375, 0.0})),
      baf::ControlPoint(std::vector<double>({1.75, 0.375, 0.0})),
      baf::ControlPoint(std::vector<double>({2.0, 0.375, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, 0.5, 0.0})),
      baf::ControlPoint(std::vector<double>({0.25, 0.5, 0.0})),
      baf::ControlPoint(std::vector<double>({0.5, 0.5, 0.0})),
      baf::ControlPoint(std::vector<double>({0.75, 0.5, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.5, 0.0})),
      baf::ControlPoint(std::vector<double>({1.25, 0.5, 0.0})),
      baf::ControlPoint(std::vector<double>({1.5, 0.5, 0.0})),
      baf::ControlPoint(std::vector<double>({1.75, 0.5, 0.0})),
      baf::ControlPoint(std::vector<double>({2.0, 0.5, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, 0.625, 0.0})),
      baf::ControlPoint(std::vector<double>({0.25, 0.625, 0.0})),
      baf::ControlPoint(std::vector<double>({0.5, 0.625, 0.0})),
      baf::ControlPoint(std::vector<double>({0.75, 0.625, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.625, 0.0})),
      baf::ControlPoint(std::vector<double>({1.25, 0.625, 0.0})),
      baf::ControlPoint(std::vector<double>({1.5, 0.625, 0.0})),
      baf::ControlPoint(std::vector<double>({1.75, 0.625, 0.0})),
      baf::ControlPoint(std::vector<double>({2.0, 0.625, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, 0.75, 0.0})),
      baf::ControlPoint(std::vector<double>({0.25, 0.75, 0.0})),
      baf::ControlPoint(std::vector<double>({0.5, 0.75, 0.0})),
      baf::ControlPoint(std::vector<double>({0.75, 0.75, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.75, 0.0})),
      baf::ControlPoint(std::vector<double>({1.25, 0.75, 0.0})),
      baf::ControlPoint(std::vector<double>({1.5, 0.75, 0.0})),
      baf::ControlPoint(std::vector<double>({1.75, 0.75, 0.0})),
      baf::ControlPoint(std::vector<double>({2.0, 0.75, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, 0.875, 0.0})),
      baf::ControlPoint(std::vector<double>({0.25, 0.875, 0.0})),
      baf::ControlPoint(std::vector<double>({0.5, 0.875, 0.0})),
      baf::ControlPoint(std::vector<double>({0.75, 0.875, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 0.875, 0.0})),
      baf::ControlPoint(std::vector<double>({1.25, 0.875, 0.0})),
      baf::ControlPoint(std::vector<double>({1.5, 0.875, 0.0})),
      baf::ControlPoint(std::vector<double>({1.75, 0.875, 0.0})),
      baf::ControlPoint(std::vector<double>({2.0, 0.875, 0.0})),

      baf::ControlPoint(std::vector<double>({0.0, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.25, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.5, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({0.75, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.25, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.5, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({1.75, 1.0, 0.0})),
      baf::ControlPoint(std::vector<double>({2.0, 1.0, 0.0})),
  };

  std::array<std::shared_ptr<baf::KnotVector>, 2> kv_ptr = {std::make_shared<baf::KnotVector>(knot_vector[0]),
                                                            std::make_shared<baf::KnotVector>(knot_vector[1])};
  std::shared_ptr<spl::NURBS<2>> nurbs_ = std::make_shared<spl::NURBS<2>>(kv_ptr, degree, control_points, weights);

  iga::LinearEquationAssembler linear_equation_assembler = iga::LinearEquationAssembler(nurbs_);
  iga::ElementIntegralCalculator elm_itg_calc = iga::ElementIntegralCalculator(nurbs_);
  int n = nurbs_->GetNumberOfControlPoints();
  std::shared_ptr<arma::dmat> matA = std::make_shared<arma::dmat>(n, n, arma::fill::zeros);
  std::shared_ptr<arma::dvec> vecB = std::make_shared<arma::dvec>(n, arma::fill::zeros);
  std::shared_ptr<arma::dvec> srcCp = std::make_shared<arma::dvec>(n, arma::fill::zeros);
  iga::itg::IntegrationRule rule = iga::itg::FourPointGaussLegendre();
};

#endif  // TEST_IGA_TEST_SPLINE_H_
