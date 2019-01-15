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
#include "vector_utils.h"

namespace spl {
template<int DIM>
class BSpline : public Spline<DIM> {
 public:
  BSpline(KnotVectors<DIM> knot_vector, std::array<Degree, DIM> degree,
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

  BSpline(const BSpline<DIM> &bspline) : Spline<DIM>(bspline) {
    PhysicalSpace<DIM> physical_space(*bspline.GetPhysicalSpace());
    physical_space_ = std::make_shared<PhysicalSpace<DIM>>(physical_space);
  }

  virtual ~BSpline() = default;

  int GetNumberOfControlPoints() override {
    return physical_space_->GetNumberOfControlPoints();
  }

  std::array<int, DIM> GetPointsPerDirection() override {
    return physical_space_->GetNumberOfPointsInEachDirection();
  }

  int GetDimension() const override {
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

  bool RemoveControlPoints(std::vector<double> scaling, int first, int last, int dimension, double tolerance) override {
    int off = first - 1, i = first, j = last;
    std::vector<double> temp = GetTemporaryNewControlPoints(scaling, first, last, off, i, j, dimension);
    if (!IsKnotRemovable(scaling[i - off - 1], temp, tolerance, i, j, off, dimension)) {
      return false;
    }
    SetNewControlPoints(temp, last, i - off, off, dimension);
    physical_space_->RemoveControlPoints(
        physical_space_->GetNumberOfControlPoints() / physical_space_->GetNumberOfPointsInEachDirection()[dimension]);
    physical_space_->DecrementNumberOfPoints(dimension);
    return true;
  }

  std::array<std::shared_ptr<spl::BSpline<DIM>>, 2> SudivideSpline(ParamCoord param_coord, int dimension) {
    this->InsertKnot(param_coord, dimension,
                     this->GetDegree(dimension).get() + 1
                         - this->GetKnotVector(dimension)->GetMultiplicity(param_coord));
    std::array<KnotVectors<DIM>, 2>
        new_knot_vectors = this->parameter_space_->GetDividedKnotVectors(param_coord, dimension);
    std::array<Degree, DIM> degrees = this->parameter_space_->GetDegrees();
    std::array<std::shared_ptr<spl::BSpline<DIM>>, 2> subdivided_splines;
    int first = 0;
    for (int i = 0; i < 2; ++i) {
      int length = new_knot_vectors[i][dimension]->GetNumberOfKnots() - degrees[dimension].get() - 1;
      std::vector<baf::ControlPoint> points = physical_space_->GetDividedControlPoints(first, length, dimension);
      spl::BSpline<DIM> spline(new_knot_vectors[i], degrees, points);
      subdivided_splines[i] = std::make_shared<spl::BSpline<DIM>>(spline);
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

  void SetNewControlPoints(std::vector<double> temp, int last, int ii, int off, int dimension) {
    std::array<int, DIM> point_handler_length = physical_space_->GetNumberOfPointsInEachDirection();
    util::MultiIndexHandler<DIM> point_handler(point_handler_length);
    std::vector<double> coordinates(GetDimension(), 0);
    int count = 0;
    for (int m = 0; m < point_handler.Get1DLength(); ++m, ++point_handler) {
      int k_index = point_handler.GetIndices()[dimension];
      int k = k_index - off;
      if (k >= 1 && k != ii && k < last - off + 2) {
        std::cout << std::endl;
        for (int l = 0; l < GetDimension(); ++l) {
          coordinates[l] = temp[(count * (last - off + 2) + k) * GetDimension() + l];
        }
        if (k == last - off + 1) count++;
        auto indices = point_handler.GetIndices();
        indices[dimension] = k < ii ? k + off : k + off - 1;
        baf::ControlPoint cp(coordinates);
        physical_space_->SetControlPoint2(indices, cp, dimension);
      }
      if ((k_index <= off && k < 1)
          || (k_index >= last + 1 && k_index < physical_space_->GetNumberOfPointsInEachDirection()[dimension])) {
        auto indices = point_handler.GetIndices();
        indices[dimension] = k_index <= off ? k_index : k_index - 1;
        physical_space_->SetControlPoint2(indices,
                                          physical_space_->GetControlPoint(point_handler.GetIndices()),
                                          dimension);
      }
    }
  }

  std::vector<double> GetTemporaryNewControlPoints(std::vector<double> scaling, int first, int last,
                                                   int &off, int &i, int &j, int dimension) const {
    std::array<int, DIM> point_handler_length = physical_space_->GetNumberOfPointsInEachDirection();
    util::MultiIndexHandler<DIM> point_handler(point_handler_length);
    int new_control_points =
        physical_space_->GetNumberOfControlPoints() / physical_space_->GetNumberOfPointsInEachDirection()[dimension];
    std::vector<double> temp(new_control_points * GetDimension() * (last - first + 3), 0);

    for (int l = 0; l < point_handler.Get1DLength(); ++l, ++point_handler) {
      if (point_handler.GetIndices()[dimension] == first - 1) {
        for (int k = 0; k < GetDimension(); ++k) {
          int index =
              point_handler.Get1DIndex() / point_handler_length[dimension] * GetDimension() * (last - first + 3);
          temp[index + k] = physical_space_->GetControlPoint(point_handler.GetIndices()).GetValue(k);
        }
      } else if (point_handler.GetIndices()[dimension] == last + 1) {
        for (int k = 0; k < GetDimension(); ++k) {
          int index =
              point_handler.Get1DIndex() / point_handler_length[dimension] * GetDimension() * (last - first + 3);
          temp[index + (last + 1 - off) * GetDimension() + k] =
              physical_space_->GetControlPoint(point_handler.GetIndices()).GetValue(k);
        }
      }
    }
    while (j - i > 0) {
      point_handler.SetIndices({0});
      for (int l = 0; l < point_handler.Get1DLength(); ++l, ++point_handler) {
        double alfi = scaling[i - first];
        double alfj = scaling[j - first];
        for (int k = 0; k < GetDimension(); ++k) {
          int index =
              point_handler.Get1DIndex() / point_handler_length[dimension] * GetDimension() * (last - first + 3);
          temp[index + (i - off) * GetDimension() + k] = (physical_space_->GetControlPoints()[i * GetDimension() + k]
              - (1 - alfi) * temp[(i - off - 1) * GetDimension() + k]) / alfi;
          temp[index + (j - off) * GetDimension() + k] =
              (physical_space_->GetControlPoints()[j * GetDimension() + k]
                  - alfj * temp[(j - off + 1) * GetDimension() + k])
                  / (1 - alfj);
        }
      }
      ++i, --j;
    }
    return temp;
  }

  bool IsKnotRemovable(double alfi, std::vector<double> temp, double tolerance,
                       int i, int j, int off, int dimension) const {
    std::array<int, DIM> point_handler_length = physical_space_->GetNumberOfPointsInEachDirection();
    point_handler_length[dimension] = 0;
    util::MultiIndexHandler<DIM> point_handler(point_handler_length);
    int new_control_points =
        physical_space_->GetNumberOfControlPoints() / physical_space_->GetNumberOfPointsInEachDirection()[dimension];
    size_t temp_length = temp.size() / new_control_points;
    for (int l = 0; l < new_control_points; ++l, ++point_handler) {
      std::vector<double> current_temp(temp_length, 0);
      for (size_t k = 0; k < temp_length; ++k) {
        current_temp[k] = temp[l * temp_length + k];
      }
      std::vector<double> temp1(GetDimension() + 1, 1);
      std::vector<double> temp2(GetDimension() + 1, 1);
      for (int k = 0; k < GetDimension(); ++k) {
        temp1[k] = current_temp[(i - off - 1) * GetDimension() + k];
        temp2[k] = current_temp[(j - off + 1) * GetDimension() + k];
      }
      if (util::VectorUtils<double>::ComputeDistance(temp1, temp2) > tolerance) {
        for (int k = 0; k < GetDimension(); ++k) {
          auto indices = point_handler.GetIndices();
          indices[dimension] = i;
          temp1[k] = physical_space_->GetControlPoint(indices).GetValue(k);
          temp2[k] = alfi * current_temp[(i - off + 1) * GetDimension() + k]
              + (1 - alfi) * current_temp[(i - off - 1) * GetDimension() + k];
        }
        if (util::VectorUtils<double>::ComputeDistance(temp1, temp2) > tolerance) {
          return false;
        }
      }
    }
    return true;
  }

  std::shared_ptr<PhysicalSpace<DIM>> physical_space_;
};
}  // namespace spl

#endif  // SRC_SPL_B_SPLINE_H_
