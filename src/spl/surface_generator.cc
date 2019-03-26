#include <utility>

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
  this->parameter_space_ = JoinParameterSpaces(nurbs_T, nurbs_C);
  this->physical_space_ = JoinPhysicalSpaces(nurbs_T, nurbs_C);
}

std::shared_ptr<spl::ParameterSpace<2>>
spl::SurfaceGenerator::JoinParameterSpaces(std::shared_ptr<spl::NURBS<1>> const &nurbs_T,
                                           std::shared_ptr<spl::NURBS<1>> const &nurbs_C) const {
  std::array<std::shared_ptr<baf::KnotVector>, 2>
      joined_knot_vector = {nurbs_T->GetKnotVector(0), nurbs_C->GetKnotVector(0)};
  std::array<Degree, 2> joined_degree = {nurbs_T->GetDegree(0), nurbs_C->GetDegree(0)};
  return std::make_shared<ParameterSpace<2>>(ParameterSpace<2>(joined_knot_vector, joined_degree));
}

std::shared_ptr<spl::WeightedPhysicalSpace<2>> spl::SurfaceGenerator::JoinPhysicalSpaces(
    std::shared_ptr<spl::NURBS<1>> const &nurbs_T,
    std::shared_ptr<spl::NURBS<1>> const &nurbs_C) const {
  std::array<int, 2> j_number_of_points =
      {nurbs_T->GetNumberOfControlPoints(), nurbs_C->GetNumberOfControlPoints()};
  std::vector<baf::ControlPoint> j_control_points;
  std::vector<double> joined_weights;
  for (int i = 0; i < nurbs_C->GetNumberOfControlPoints(); ++i) {
    std::array<int, 1> index_space_2 = {i};
    for (int j = 0; j < nurbs_T->GetNumberOfControlPoints(); ++j) {
      std::array<int, 1> index_space_1 = {j};
      j_control_points.emplace_back(nurbs_T->GetControlPoint(index_space_1) +
          nurbs_C->GetControlPoint(index_space_2));
      joined_weights.emplace_back(nurbs_T->GetWeight(index_space_1) * nurbs_C->GetWeight(index_space_2));
    }
  }
  return std::make_shared<spl::WeightedPhysicalSpace<2>>(spl::WeightedPhysicalSpace<2>(
      j_control_points, joined_weights, j_number_of_points));
}

spl::SurfaceGenerator::SurfaceGenerator(std::shared_ptr<spl::NURBS<1>> const &nurbs_T,
                                        std::shared_ptr<spl::NURBS<1>> const &nurbs_C,
                                        int nbInter, std::vector<std::array<double, 3>> scaling) {
  double step_size = nurbs_T->GetKnotVectorRange(0) / (nbInter - 1);
  std::vector<ParamCoord> v_i;
  std::vector<double> dT_v;
  std::vector<double> ddT_v;
  std::vector<double> t_v;
  const std::vector<int> dimensions = {0, 1, 2};
  std::array<int, 1> first_derivative = {1};
  std::array<int, 1> second_derivative = {2};
  std::array<std::array<double, 4>, 4> transMatrix({std::array<double, 4>({0.0, 0.0, 0.0, 0.0}),
                                                    std::array<double, 4>({0.0, 0.0, 0.0, 0.0}),
                                                    std::array<double, 4>({0.0, 0.0, 1.0, 0.0}),
                                                    std::array<double, 4>({0.0, 0.0, 0.0, 1.0})});
  double weight_v;
  int m = nurbs_C->GetNumberOfControlPoints();
  std::vector<baf::ControlPoint> j_control_points(nbInter * m, baf::ControlPoint({0.0, 0.0, 0.0}));
  std::vector<double> j_weights(nbInter * m, 0.0);
  v_i.reserve(nbInter);
  for (int i = 0; i < nbInter; ++i) {
    v_i.emplace_back(ParamCoord{i * step_size});
    t_v = nurbs_T->Evaluate(std::array<ParamCoord, 1>({v_i[i]}), dimensions);
    weight_v = nurbs_T->Evaluate(std::array<ParamCoord, 1>({v_i[i]}), std::vector<int>({3}))[0];
    dT_v = nurbs_T->EvaluateDerivative(std::array<ParamCoord, 1>({v_i[i]}), dimensions, first_derivative);
    ddT_v = nurbs_T->EvaluateDerivative(std::array<ParamCoord, 1>({v_i[i]}), dimensions, second_derivative);
    transMatrix = GetTransformation(t_v,
                                    dT_v,
                                    ddT_v,
                                    std::array<double, 3>({transMatrix[0][2], transMatrix[1][2], transMatrix[2][2]}),
                                    i);
    for (int j = 0; j < m; ++j) {
      int indexControlPoint = j * nbInter + i;
      baf::ControlPoint control_point_j = nurbs_C->GetControlPoint(std::array<int, 1>({j}));
      j_control_points[indexControlPoint] = control_point_j.Transform(transMatrix, scaling[i]);
      j_weights[indexControlPoint] = nurbs_C->GetWeight(std::array<int, 1>({j})) * weight_v;
    }
  }
  baf::KnotVector knot_vector_t(v_i, nurbs_T->GetDegree(0), nbInter);
  std::shared_ptr<baf::KnotVector> knot_vector_t_ptr = std::make_shared<baf::KnotVector>(knot_vector_t);
  std::array<std::shared_ptr<baf::KnotVector>, 2> joined_knot_vector =
      {knot_vector_t_ptr, nurbs_C->GetKnotVector(0)};
  std::array<int, 2> j_number_of_points = {nbInter, m};
  std::array<Degree, 2> joined_degree = {nurbs_T->GetDegree(0), nurbs_C->GetDegree(0)};
  this->parameter_space_ = std::make_shared<ParameterSpace<2>>(ParameterSpace<2>(
      joined_knot_vector, joined_degree));
  this->physical_space_ = std::make_shared<spl::WeightedPhysicalSpace<2>>(spl::WeightedPhysicalSpace<2>(
      j_control_points, j_weights, j_number_of_points));
}

std::array<double, 3> spl::SurfaceGenerator::CrossProduct(std::vector<double> a, std::vector<double> b) const {
  std::array<double, 3> result = {a[1] * b[2] - a[2] * b[1],
                                  a[2] * b[0] - a[0] * b[2],
                                  a[0] * b[1] - a[1] * b[0]};
  return result;
}

std::array<double, 3> spl::SurfaceGenerator::CrossProduct(std::array<double, 3> a, std::array<double, 3> b) const {
  std::array<double, 3> result = {a[1] * b[2] - a[2] * b[1],
                                  a[2] * b[0] - a[0] * b[2],
                                  a[0] * b[1] - a[1] * b[0]};
  return result;
}

double spl::SurfaceGenerator::DotProduct(std::array<double, 3> a, std::array<double, 3> b) const {
  double dotProduct = 0.0;
  for (int i = 0; i < 3; ++i) {
    dotProduct += a[i] * b[i];
  }
  return dotProduct;
}

double spl::SurfaceGenerator::DotProduct(std::array<double, 3> a, std::vector<double> b) const {
  double dotProduct = 0.0;
  for (int i = 0; i < 3; ++i) {
    dotProduct += a[i] * b[i];
  }
  return dotProduct;
}

double spl::SurfaceGenerator::ComputeNorm(std::vector<double> a) {
  double norm = 0.0;
  for (double i : a) {
    norm += std::pow(i, 2);
  }
  return std::sqrt(norm);
}

double spl::SurfaceGenerator::ComputeNorm(std::array<double, 3> a) {
  double norm = 0.0;
  for (double i : a) {
    norm += std::pow(i, 2);
  }
  return std::sqrt(norm);
}

std::array<double, 3> spl::SurfaceGenerator::ComputeNormal(std::vector<double> T,
                                                           std::vector<double> dT,
                                                           std::vector<double> ddT,
                                                           std::array<double, 3> previous,
                                                           int index) {
  std::array<double, 3> crossProduct = CrossProduct(dT, std::move(ddT));
  double norm_cross = ComputeNorm(crossProduct);
  if (norm_cross > 0.0) {
    for (int i = 0; i < 3; ++i) {
      crossProduct[i] /= norm_cross;
    }
    if (DotProduct(crossProduct, previous) == 0.0 && index != 0) {
      throw std::runtime_error(
          "The swept surface is discontinous due to a rotation. "
          "Consider increasing the degree of the trajectory NURBS.");
    }
    return crossProduct;
  }
  double norm_dT = ComputeNorm(dT);
  for (int i = 0; i < 3; ++i) {
    dT[i] /= norm_dT;
  }
  double product = DotProduct(previous, T);
  for (int i = 0; i < 3; ++i) {
    crossProduct[i] = previous[i] - product * T[i];
  }
  norm_cross = ComputeNorm(crossProduct);
  for (int i = 0; i < 3; ++i) {
    crossProduct[i] /= norm_cross;
  }
  return crossProduct;
}

std::array<std::array<double, 4>, 4> spl::SurfaceGenerator::GetTransformation(std::vector<double> t,
                                                                              std::vector<double> dT,
                                                                              std::vector<double> ddT,
                                                                              std::array<double, 3> prev_z,
                                                                              int index) {
  std::array<std::array<double, 4>, 4> transformation({std::array<double, 4>({0.0, 0.0, 0.0, 0.0}),
                                                       std::array<double, 4>({0.0, 0.0, 0.0, 0.0}),
                                                       std::array<double, 4>({0.0, 0.0, 0.0, 0.0}),
                                                       std::array<double, 4>({0.0, 0.0, 0.0, 1.0})});
  std::array<double, 3> z = ComputeNormal(t, dT, std::move(ddT), prev_z, index);
  double norm_dT = 0.0;
  for (int j = 0; j < 3; ++j) {
    norm_dT += std::pow(dT[j], 2);
    transformation[j][3] = t[j];
  }
  norm_dT = std::sqrt(norm_dT);
  for (int j = 0; j < 3; ++j) {
    transformation[j][0] = dT[j] / norm_dT;
    transformation[j][2] = z[j];
  }
  std::array<double, 3> x({transformation[0][0], transformation[1][0], transformation[2][0]});
  std::array<double, 3> y = CrossProduct(z, x);
  for (int i = 0; i < 3; ++i) {
    transformation[i][1] = y[i];
  }
  return transformation;
}



