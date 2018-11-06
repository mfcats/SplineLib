/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "surface_generator.h"

spl::SurfaceGenerator::SurfaceGenerator(std::shared_ptr<spl::NURBS<1>> const &nurbs_T,
                                        std::shared_ptr<spl::NURBS<1>> const &nurbs_C) {
  this->parameter_space_ = JoinParameterSpaces(nurbs_T->GetParameterSpace(), nurbs_C->GetParameterSpace());
  this->physical_space_ = JoinPhysicalSpaces(nurbs_T->GetPhysicalSpace(), nurbs_C->GetPhysicalSpace());
}

std::shared_ptr<spl::ParameterSpace<2>> spl::SurfaceGenerator::JoinParameterSpaces(
    std::shared_ptr<ParameterSpace<1>> const &space_1, std::shared_ptr<ParameterSpace<1>> const &space_2) const {
  std::array<std::shared_ptr<baf::KnotVector>, 2> joined_knot_vector =
      {space_1->GetKnotVector(0), space_2->GetKnotVector(0)};
  std::array<Degree, 2> joined_degree = {space_1->GetDegree(0), space_2->GetDegree(0)};
  return std::make_shared<ParameterSpace<2>>(ParameterSpace<2>(joined_knot_vector, joined_degree));
}

std::shared_ptr<spl::WeightedPhysicalSpace<2>> spl::SurfaceGenerator::JoinPhysicalSpaces(
    std::shared_ptr<spl::PhysicalSpace<1>> const &space_1,
    std::shared_ptr<spl::PhysicalSpace<1>> const &space_2) const {
  std::array<int, 2> j_number_of_points =
      {space_1->GetNumberOfControlPoints(), space_2->GetNumberOfControlPoints()};
  std::vector<baf::ControlPoint> j_control_points;
  std::vector<double> joined_weights;
  for (int i = 0; i < space_2->GetNumberOfControlPoints(); ++i) {
    std::array<int, 1> index_space_2 = {i};
    for (int j = 0; j < space_1->GetNumberOfControlPoints(); ++j) {
      std::array<int, 1> index_space_1 = {j};
      j_control_points.emplace_back(space_1->GetControlPoint(index_space_1) +
                                         space_2->GetControlPoint(index_space_2));
      joined_weights.emplace_back(space_1->GetWeight(index_space_1) * space_2->GetWeight(index_space_2));
    }
  }
  return std::make_shared<spl::WeightedPhysicalSpace<2>>(spl::WeightedPhysicalSpace<2>(
      j_control_points, joined_weights, j_number_of_points));
}

spl::SurfaceGenerator::SurfaceGenerator(std::shared_ptr<spl::NURBS<1>> const &nurbs_T,
                                        std::shared_ptr<spl::NURBS<1>> const &nurbs_C,
                                        int nbInter, std::vector<double> scaling) {
  // ToDo: stepsize, weighting function
  double step_size = nurbs_T->GetKnotVectorRange(1)/(nbInter-1);
  std::vector<ParamCoord> v_i;
  std::vector<double> dT_v;
  std::vector<double> ddT_v;
  std::vector<double> t_v;
  std::array<double, 3> cross_product;
  const std::vector<int> dimensions = {0, 1, 2};
  std::array<int, 1> first_derivative = {1};
  std::array<int, 1> second_derivative = {2};

  std::array<double, 3> b_v;
  std::array<double, 3> x_v;
  std::array<double, 3> y_v;
  std::array<double, 3> z_v;
  std::array<double, 3> o_v;

  double norm_cross = 0.0;
  double norm_dT_v = 0.0;

  std::vector<baf::ControlPoint> j_control_points;
  std::vector<double> joined_weights;
  joined_weights.reserve(nbInter);

  for (int i = 0; i < nbInter; ++i) {
    v_i.emplace_back(ParamCoord{i * step_size});
    t_v = nurbs_T->Evaluate(v_i[i], dimensions);
    dT_v = nurbs_T->EvaluateDerivative(v_i[i], dimensions, first_derivative);
    ddT_v = nurbs_T->EvaluateDerivative(v_i[i], dimensions, second_derivative);
    cross_product = CrossProduct(dT_v, ddT_v):
    for (int j = 0; j < 3; ++j) {
      norm_cross += std::pow(cross_product[j], 2);
      norm_dT_v += std::pow(dT_v[j], 2);
    }
    norm_cross = std::sqrt(norm_cross);
    norm_dT_v = std::sqrt(norm_dT_v);
    for (int j = 0; j < 3; ++j) {
      b_v[j] = cross_product[j]/norm_cross;
      x_v[j] = dT_v[j]/norm_dT_v;
      z_v[j] = b_v[j];
    }
    y_v = CrossProduct(z_v, x_v);
    for (int j = 0; j < m; ++j) {
      baf::ControlPoint control_point_j = nurbs_C->GetControlPoint(j);
      j_control_points.emplace_back(control_point_j.Transform(x_v, y_v, z_v, o_v));
    }
    joined_weights.push_back(InterpolateWeights(nurbs_T->GetWeights(), v_i[i]));
  }


  for (int i = 0; i < m; ++i) {
    std::array<int, 1> index_space_2 = {i};
    for (int j = 0; j < n; ++j) {
      std::array<int, 1> index_space_1 = {j};
      joined_weights.emplace_back(nurbs_T->GetPhysicalSpace()->GetWeight(index_space_1) *
      nurbs_C->GetPhysicalSpace()->GetWeight(index_space_2));
      j_control_points.emplace_back(nurbs_T->GetPhysicalSpace()->GetControlPoint(index_space_1) +
          nurbs_C->GetPhysicalSpace()->GetControlPoint(index_space_2));
    }
  }
  // Divide by weights
  std::shared_ptr<baf::KnotVector> knot_vector_t_ptr = std::make_shared<baf::KnotVector>(knot_vector_t);
  std::array<std::shared_ptr<baf::KnotVector>, 2> joined_knot_vector =
      {knot_vector_t_ptr, nurbs_C->GetParameterSpace()->GetKnotVector(0)};
  std::array<int, 2> j_number_of_points = {m, n};
  std::array<Degree, 2> joined_degree = {nurbs_T->GetParameterSpace()->GetDegree(0), nurbs_C->GetParameterSpace()->GetDegree(0)};
  this->parameter_space_ = std::make_shared<ParameterSpace<2>>(ParameterSpace<2>(
      joined_knot_vector, joined_degree));
  this->physical_space_ = std::make_shared<spl::WeightedPhysicalSpace<2>>(spl::WeightedPhysicalSpace<2>(
      j_control_points, joined_weights, j_number_of_points));
}

double spl::SurfaceGenerator::ComputeNormal(
    std::shared_ptr<spl::NURBS<1>> const &nurbs, ParamCoord param_coord, int direction) const {
  // Evaluates the formula 10.25 in the NURBS-book
  double normal_value = 1.0;
  return normal_value;
}

std::array<double, 3> spl::SurfaceGenerator::CrossProduct( std::vector<double> a, std::vector<double> b) const {
  std::array<double, 3> result = {a[1] * b[2] - a[2] * b[1],
                                  a[2] * b[0] - a[0] * b[2],
                                  a[0] * b[1] - a[1] * b[0]};
  return result;
}

std::array<double, 3> spl::SurfaceGenerator::CrossProduct( std::array<double, 3> a, std::array<double, 3> b) const {
  std::array<double, 3> result = {a[1] * b[2] - a[2] * b[1],
                                  a[2] * b[0] - a[0] * b[2],
                                  a[0] * b[1] - a[1] * b[0]};
  return result;
}

baf::KnotVector spl::SurfaceGenerator::AverageKnots(Degree degree, int nbControlPoints) {
  std::vector<ParamCoord> knots;
  for (int i = 0; i < degree.get(); ++i) {
    knots.emplace_back(ParamCoord{0.0});
  }
  double tempParam;
  for (int i = 1; i <= nbControlPoints - 1 - degree.get(); ++i) {
    tempParam = 0;
    for (int j = i; j < i + degree.get(); ++j) {
      tempParam += u[j];
    }
    tempParam /= degree.get();
    knots.emplace_back(ParamCoord{tempParam});
  }
  for (int i = 0; i < degree.get(); ++i) {
    knots.emplace_back(ParamCoord{1.0});
  }
  return baf::KnotVector(knots);
}

double spl::SurfaceGenerator::InterpolateWeights(std::shared_ptr<spl::NURBS<1>> const &nurbs,
    std::vector<ParamCoord> param_coords) {
  

}


