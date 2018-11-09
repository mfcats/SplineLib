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
#include <iostream>
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
    std::array<int, DIM> number_of_points = physical_space_->GetNumberOfPointsInEachDirection();
    for (auto &number : number_of_points) {
      --number;
    }
    util::MultiIndexHandler<DIM> point_handler(physical_space_->GetNumberOfPointsInEachDirection());
    point_handler.SetIndices(number_of_points);
    for (int i = point_handler.Get1DLength() - 1; i >= 0; --i, --point_handler) {
      auto current_point = point_handler.GetIndices()[dimension];
      if (current_point <= last && point_handler.GetIndices()[dimension] >= first) {
        std::array<int, DIM> indices0 = point_handler.GetIndices(), indices1 = indices0;
        --indices1[dimension];
        baf::ControlPoint cp0 = physical_space_->GetControlPoint(indices0);
        baf::ControlPoint cp1 = physical_space_->GetControlPoint(indices1);
        std::vector<double> coordinates;
        for (int j = 0; j < cp0.GetDimension(); ++j) {
          coordinates.push_back(scaling[current_point - first] * cp0.GetValue(j)
                                    + (1 - scaling[current_point - first]) * cp1.GetValue(j));
        }
        baf::ControlPoint new_cp(coordinates);
        current_point != last ? physical_space_->SetControlPoint(indices0, new_cp)
                              : physical_space_->InsertControlPoint(indices0, new_cp);
      }
    }
    physical_space_->IncrementNumberOfPoints(dimension);
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

  std::shared_ptr<PhysicalSpace<DIM>> physical_space_;
};
}  // namespace spl

#endif  // SRC_SPL_B_SPLINE_H_
