/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_B_SPLINE_H_
#define SRC_SPL_B_SPLINE_H_

#include <algorithm>
#include <array>
#include <functional>
#include <vector>

#include "b_spline_generator.h"
#include "spline.h"
#include "spline_generator.h"

namespace spl {
template<int DIM>
class BSpline : public Spline<DIM> {
 public:
  BSpline(std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vector,
          std::array<Degree, DIM> degree,
          const std::vector<baf::ControlPoint> &control_points) : Spline<DIM>(knot_vector, degree) {
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = knot_vector[i]->GetNumberOfKnots() - degree[i].get() - 1;
    }
    physical_space_ = std::make_shared<PhysicalSpace<DIM>>(PhysicalSpace<DIM>(control_points, number_of_points));
  }

  explicit BSpline(BSplineGenerator<DIM> b_spline_generator) : Spline<DIM>(b_spline_generator.GetParameterSpace()) {
    physical_space_ = b_spline_generator.GetPhysicalSpace();
  }

  virtual ~BSpline() = default;

  BSpline(ParameterSpace<DIM> parameter_space, PhysicalSpace<DIM> physical_space)
      : Spline<DIM>(std::make_shared<spl::ParameterSpace<DIM>>(parameter_space)),
        physical_space_(std::make_shared<spl::PhysicalSpace<DIM>>(physical_space)) {
  }

  int GetNumberOfControlPoints() override {
    return physical_space_->GetNumberOfControlPoints();
  }

  std::array<int, DIM> GetPointsPerDirection() override {
    return physical_space_->GetNumberOfPointsInEachDirection();
  }

  int GetDimension() override {
    return physical_space_->GetDimension();
  }

  double GetControlPoint(std::array<int, DIM> indices, int dimension) override {
    return physical_space_->GetControlPoint(indices).GetValue(dimension);
  }

  std::shared_ptr<spl::PhysicalSpace<DIM>> GetPhysicalSpace() const override {
    return physical_space_;
  }

  void AdjustControlPoints(std::vector<double> scaling, int first, int last, int dimension) override {
    std::array<int, DIM> point_handler_length = physical_space_->GetNumberOfPointsInEachDirection();
    ++point_handler_length[dimension];
    util::MultiIndexHandler<DIM> point_handler(point_handler_length);
    std::array<int, DIM> maximum_point_index = physical_space_->GetMaximumPointIndexInEachDirection();
    ++maximum_point_index[dimension];
    point_handler.SetIndices(maximum_point_index);
    physical_space_->AddControlPoints(physical_space_->GetNumberOfControlPoints() / maximum_point_index[dimension]);
    for (int i = point_handler.Get1DLength() - 1; i >= 0; --i, --point_handler) {
      auto current_point = point_handler.GetIndices()[dimension];
      std::array<int, DIM> indices = point_handler.GetIndices();
      baf::ControlPoint new_control_point = GetNewControlPoint(indices, dimension, scaling, current_point, first, last);
      physical_space_->SetControlPoint(indices, new_control_point, dimension);
    }
    physical_space_->IncrementNumberOfPoints(dimension);
  }

  std::array<std::any, 2>
  GetSubdividedSpline(std::array<std::array<std::shared_ptr<baf::KnotVector>, DIM>, 2> knot_vectors,
                      int dimension, std::array<Degree, DIM> degrees) override {
    std::array<std::any, 2> subdivided_splines;
    int first = 0;
    for (int i = 0; i < 2; ++i) {
      int length = knot_vectors[i][dimension]->GetNumberOfKnots() - degrees[dimension].get() - 1;
      std::vector<baf::ControlPoint> points = physical_space_->GetSplittedControlPoints(first, length, dimension);
      spl::BSpline<DIM> spline(knot_vectors[i], degrees, points);
      std::shared_ptr<spl::BSpline<DIM>> spline_ptr = std::make_shared<spl::BSpline<DIM>>(spline);
      subdivided_splines[i] = std::make_any<std::shared_ptr<spl::BSpline<DIM>>>(spline_ptr);
      first = length;
    }
    return subdivided_splines;
  }

 private:
  double GetEvaluatedControlPoint(std::array<ParamCoord, DIM> param_coord,
                                  std::array<int, DIM> indices,
                                  int dimension) const override {
    return this->parameter_space_->GetBasisFunctions(indices, param_coord)
        * physical_space_->GetControlPoint(indices).GetValue(dimension);
  }

  double GetEvaluatedDerivativeControlPoint(std::array<ParamCoord, DIM> param_coord,
                                            std::array<int, DIM> derivative,
                                            std::array<int, DIM> indices,
                                            int dimension) const override {
    return this->parameter_space_->GetBasisFunctionDerivatives(indices, param_coord, derivative)
        * physical_space_->GetControlPoint(indices).GetValue(dimension);
  }

  baf::ControlPoint GetNewControlPoint(std::array<int, DIM> indices, int dimension, std::vector<double> scaling,
                                       int current_point, int first, int last) {
    if (current_point > last) {
      --indices[dimension];
      return physical_space_->GetControlPoint(indices);
    } else if (current_point >= first) {
      std::array<int, DIM> lower_indices = indices;
      --lower_indices[dimension];
      baf::ControlPoint upper_control_point = physical_space_->GetControlPoint(indices);
      baf::ControlPoint lower_control_point = physical_space_->GetControlPoint(lower_indices);
      std::vector<double> coordinates;
      for (int j = 0; j < upper_control_point.GetDimension(); ++j) {
        coordinates.push_back(scaling[current_point - first] * upper_control_point.GetValue(j)
                                  + (1 - scaling[current_point - first]) * lower_control_point.GetValue(j));
      }
      return baf::ControlPoint(coordinates);
    } else {
      return physical_space_->GetControlPoint(indices);
    }
  }

  std::shared_ptr<PhysicalSpace<DIM>> physical_space_;
};
}  // namespace spl

#endif  // SRC_SPL_B_SPLINE_H_
