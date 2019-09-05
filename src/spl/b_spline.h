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
#include <utility>
#include <vector>

#include "b_spline_generator.h"
#include "multi_index_handler.h"
#include "numeric_operations.h"
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

  bool AreEqual(const BSpline<DIM> &rhs, double tolerance = util::NumericSettings<double>::kEpsilon()) const {
    return this->parameter_space_->AreEqual(*rhs.parameter_space_.get(), tolerance)
        && physical_space_->AreEqual(*rhs.physical_space_.get(), tolerance);
  }

  void AdjustControlPoints(std::vector<double> scaling, int first, int last, int dimension) override {
    std::array<int, DIM> point_handler_length = this->GetPointsPerDirection();
    ++point_handler_length[dimension];
    util::MultiIndexHandler<DIM> point_handler(point_handler_length);
    std::array<int, DIM> maximum_point_index = physical_space_->GetMaximumPointIndexInEachDirection();
    ++maximum_point_index[dimension];
    point_handler.SetIndices(maximum_point_index);
    physical_space_->AddControlPoints(this->GetNumberOfControlPoints() / maximum_point_index[dimension]);
    for (int i = point_handler.Get1DLength() - 1; i >= 0; --i, --point_handler) {
      auto current_point_index = point_handler[dimension];
      std::array<int, DIM> indices = point_handler.GetIndices();
      baf::ControlPoint new_control_point = GetNewControlPoint(indices, dimension, scaling,
                                                               current_point_index, first, last);
      physical_space_->SetControlPoint(indices, new_control_point, dimension,
                                       util::NumericOperations<int>::increment);
    }
    physical_space_->IncrementNumberOfPoints(dimension);
  }

  bool RemoveControlPoints(std::vector<double> scaling, int first, int last, int dimension, double tolerance) override {
    int off = first - 1, i = first, j = last;
    std::vector<double> temp = GetTempNewControlPoints(scaling, off, last, i, j, dimension);
    auto diff = static_cast<int>(ceil((j - i) / 2.0));
    i += diff, j -= diff;
    if (!IsKnotRemovable(scaling[i - off - 1], temp, tolerance, i, j, off, dimension)) {
      return false;
    }
    SetNewControlPoints(temp, last, i - off, off, dimension);
    physical_space_->RemoveControlPoints(this->GetNumberOfControlPoints() / this->GetPointsPerDirection()[dimension]);
    physical_space_->DecrementNumberOfPoints(dimension);
    return true;
  }

  std::array<std::shared_ptr<spl::BSpline<DIM>>, 2> SudivideSpline(ParamCoord param_coord, int dimension) {
    this->InsertKnot(param_coord, dimension, this->GetDegree(dimension).get() + 1
                         - this->GetKnotVector(dimension)->GetMultiplicity(param_coord));
    std::array<KnotVectors<DIM>, 2> new_knots = this->parameter_space_->GetDividedKnotVectors(param_coord, dimension);
    std::array<Degree, DIM> degrees;
    for (int i = 0; i < DIM; ++i) {
      degrees[i] = this->GetDegree(i);
    }
    std::array<std::shared_ptr<spl::BSpline<DIM>>, 2> subdivided_splines;
    int first = 0;
    for (int i = 0; i < 2; ++i) {
      int length = new_knots[i][dimension]->GetNumberOfKnots() - this->GetDegree(dimension).get() - 1;
      std::vector<baf::ControlPoint> points = physical_space_->GetDividedControlPoints(first, length, dimension);
      spl::BSpline<DIM> spline(new_knots[i], degrees, points);
      subdivided_splines[i] = std::make_shared<spl::BSpline<DIM>>(spline);
      first = length;
    }
    return subdivided_splines;
  }

 private:
  std::shared_ptr<spl::PhysicalSpace<DIM>> GetPhysicalSpace() const override {
    return physical_space_;
  }

  double GetEvaluatedControlPoint(std::array<ParamCoord, DIM> param_coord,
                                  std::array<int, DIM> indices,
                                  int dimension) const override {
    return this->parameter_space_->GetBasisFunctions(indices, param_coord) * this->GetControlPoint(indices, dimension);
  }

  double GetEvaluatedDerivativeControlPoint(std::array<ParamCoord, DIM> param_coord,
                                            std::array<int, DIM> derivative,
                                            std::array<int, DIM> indices,
                                            int dimension) const override {
    return this->parameter_space_->GetBasisFunctionDerivatives(indices, param_coord, derivative)
        * this->GetControlPoint(indices, dimension);
  }

  baf::ControlPoint GetNewControlPoint(std::array<int, DIM> indices, int dimension, std::vector<double> scaling,
                                       int current_point_index, int first, int last) {
    if (current_point_index > last) {
      --indices[dimension];
      return this->GetControlPoint(indices);
    } else if (current_point_index >= first) {
      std::array<int, DIM> lower_indices = indices;
      --lower_indices[dimension];
      baf::ControlPoint upper_control_point = this->GetControlPoint(indices);
      baf::ControlPoint lower_control_point = this->GetControlPoint(lower_indices);
      std::vector<double> coordinates;
      for (int j = 0; j < upper_control_point.GetDimension(); ++j) {
        coordinates.push_back(scaling[current_point_index - first] * upper_control_point.GetValue(j)
                                  + (1 - scaling[current_point_index - first]) * lower_control_point.GetValue(j));
      }
      return baf::ControlPoint(coordinates);
    } else {
      return this->GetControlPoint(indices);
    }
  }

  void SetNewControlPoints(const std::vector<double> &temp, int last, int ii, int off, int dimension) {
    std::array<int, DIM> point_handler_length = this->GetPointsPerDirection();
    util::MultiIndexHandler<DIM> point_handler(point_handler_length);
    for (int m = 0; m < point_handler.Get1DLength(); ++m, ++point_handler) {
      int k = point_handler[dimension];
      if (k - off >= 1 && k - off != ii && k < last + 2) {
        int index = (point_handler.ExtractDimension(dimension) * (last - off + 2) + k - off) * this->GetPointDim();
        std::vector<double> coordinates(temp.begin() + index, temp.begin() + index + this->GetPointDim());
        auto indices = point_handler.GetIndices();
        indices[dimension] = k - off < ii ? k : k - 1;
        physical_space_->SetControlPoint(indices, baf::ControlPoint(coordinates), dimension,
                                         util::NumericOperations<int>::decrement);
      }
      if ((k <= off && k - off < 1) || (k >= last + 1 && k < this->GetPointsPerDirection()[dimension])) {
        auto indices = point_handler.GetIndices();
        indices[dimension] = k <= off ? k : k - 1;
        physical_space_->SetControlPoint(indices, this->GetControlPoint(point_handler.GetIndices()), dimension,
                                         util::NumericOperations<int>::decrement);
      }
    }
  }

  std::vector<double> GetTempNewControlPoints(const std::vector<double> &scaling, int off, int last,
                                              int i, int j, int dimension) const {
    std::array<int, DIM> point_handler_length = this->GetPointsPerDirection();
    util::MultiIndexHandler<DIM> point_handler(point_handler_length);
    int new_points = this->GetNumberOfControlPoints() / this->GetPointsPerDirection()[dimension];
    std::vector<double> temp(new_points * this->GetPointDim() * (last - off + 2), 0);
    std::shared_ptr<std::vector<double>> temp_ptr = std::make_shared<std::vector<double>>(temp);
    for (int l = 0; l < point_handler.Get1DLength(); ++l, ++point_handler) {
      if (point_handler[dimension] == off || point_handler[dimension] == last + 1) {
        int diff = point_handler[dimension] == off ? off : last + 1;
        SetTempNewControlPoint(point_handler, temp_ptr, 1, diff, off, last, dimension, 0);
      }
    }
    for (; j - i > 0; ++i, --j) {
      point_handler.SetIndices({0});
      for (int l = 0; l < point_handler.Get1DLength(); ++l, ++point_handler) {
        if (point_handler[dimension] == i) {
          SetTempNewControlPoint(point_handler, temp_ptr, scaling[i - off - 1], i, off, last, dimension, -1);
        }
        if (point_handler[dimension] == j) {
          SetTempNewControlPoint(point_handler, temp_ptr, 1 - scaling[j - off - 1], j, off, last, dimension, 1);
        }
      }
    }
    return *temp_ptr;
  }

  void SetTempNewControlPoint(const util::MultiIndexHandler<DIM> &point_handler,
                              const std::shared_ptr<std::vector<double>> &temp_ptr,
                              double alpha, int x, int off, int last, int dimension, int shift) const {
    int index = point_handler.ExtractDimension(dimension) * (last - off + 2);
    for (int k = 0; k < this->GetPointDim(); ++k) {
      (*temp_ptr)[(index + x - off) * this->GetPointDim() + k] = (this->GetControlPoint(point_handler.GetIndices(), k)
          - (1 - alpha) * (*temp_ptr)[(index + x - off + shift) * this->GetPointDim() + k]) / alpha;
    }
  }

  bool IsKnotRemovable(double alfi, const std::vector<double> &temp, double tolerance,
                       int i, int j, int off, int dimension) const {
    std::array<int, DIM> point_handler_length = this->GetPointsPerDirection();
    point_handler_length[dimension] = 0;
    util::MultiIndexHandler<DIM> point_handler(point_handler_length);
    int new_points = this->GetNumberOfControlPoints() / this->GetPointsPerDirection()[dimension];
    size_t temp_length = temp.size() / new_points;
    for (int l = 0; l < new_points; ++l, ++point_handler) {
      std::vector<double> temp1(this->GetPointDim() + 1, 1), temp2(this->GetPointDim() + 1, 1);
      for (int k = 0; k < this->GetPointDim(); ++k) {
        temp1[k] = temp[l * temp_length + (i - off - 1) * this->GetPointDim() + k];
        temp2[k] = temp[l * temp_length + (j - off + 1) * this->GetPointDim() + k];
      }
      if (util::VectorUtils<double>::ComputeDistance(temp1, temp2) > tolerance) {
        for (int k = 0; k < this->GetPointDim(); ++k) {
          auto indices = point_handler.GetIndices();
          indices[dimension] = i;
          temp1[k] = this->GetControlPoint(indices, k);
          temp2[k] = alfi * temp[l * temp_length + (i - off + 1) * this->GetPointDim() + k]
              + (1 - alfi) * temp[l * temp_length + (i - off - 1) * this->GetPointDim() + k];
        }
        if (util::VectorUtils<double>::ComputeDistance(temp1, temp2) > tolerance) {
          return false;
        }
      }
    }
    return true;
  }

  void SetNewControlPoint(baf::ControlPoint control_point, double /*weight*/, std::array<int, DIM> indices) override {
    physical_space_->SetControlPoint(indices, control_point);
  }

  std::shared_ptr<PhysicalSpace<DIM>> physical_space_;
};
}  // namespace spl

#endif  // SRC_SPL_B_SPLINE_H_
