/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_SPL_NURBS_H_
#define SRC_SPL_NURBS_H_

#include <algorithm>
#include <array>
#include <functional>
#include <utility>
#include <vector>

#include "src/spl/b_spline.h"
#include "src/spl/spline.h"
#include "src/spl/weighted_physical_space.h"

namespace splinelib::src::spl {
template<int PARAMETRIC_DIMENSIONALITY>
class NURBS : public Spline<PARAMETRIC_DIMENSIONALITY> {
 public:
  NURBS(baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> knot_vector, std::array<Degree, PARAMETRIC_DIMENSIONALITY> degree,
        const std::vector<spl::ControlPoint> &control_points,
        std::vector<Weight> weights) : Spline<PARAMETRIC_DIMENSIONALITY>(knot_vector, degree) {
    std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      number_of_points[i] = knot_vector[i]->GetNumberOfKnots() - degree[i].Get() - 1;
    }
    physical_space_ = std::make_shared<WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>>(
        WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>(control_points, weights, number_of_points));
  }

  NURBS(std::shared_ptr<WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>> weighted_physical_space,
        std::shared_ptr<ParameterSpace<PARAMETRIC_DIMENSIONALITY>> parameter_space) {
    physical_space_ = weighted_physical_space;
    this->parameter_space_ = parameter_space;
  }

  NURBS(const NURBS<PARAMETRIC_DIMENSIONALITY> &nurbs) : Spline<PARAMETRIC_DIMENSIONALITY>(nurbs) {
    WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY> weighted_physical_space(*nurbs.physical_space_);
    physical_space_ = std::make_shared<WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>>(weighted_physical_space);
  }

  NURBS(const NURBS<PARAMETRIC_DIMENSIONALITY> &nurbs, const std::vector<spl::ControlPoint> &control_points)
      : Spline<PARAMETRIC_DIMENSIONALITY>(nurbs) {
    WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>
        weighted_physical_space(control_points, nurbs.physical_space_->GetWeights(),
                                nurbs.GetNumberOfPointsPerDirection());
    physical_space_ = std::make_shared<WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>>(weighted_physical_space);
  }

  NURBS(NURBS<PARAMETRIC_DIMENSIONALITY> &&other) = delete;
  NURBS & operator=(const NURBS<PARAMETRIC_DIMENSIONALITY> &rhs) = delete;
  NURBS & operator=(NURBS<PARAMETRIC_DIMENSIONALITY> &&rhs) = delete;
  ~NURBS() override = default;

  bool AreEqual(const NURBS<PARAMETRIC_DIMENSIONALITY> &rhs, double tolerance =
      util::numeric_settings::GetEpsilon<double>()) const {
    return this->parameter_space_->AreEqual(*rhs.parameter_space_.get(), tolerance)
        && physical_space_->AreEqual(*rhs.physical_space_.get(), Tolerance{tolerance});
  }

  double GetHomogeneousControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> indices, int dimension) const {
    return physical_space_->GetHomogeneousControlPoint(indices)[Dimension{dimension}];
  }

  void AdjustControlPoints(std::vector<double> scaling, int first, int last, int dimension) override {
    physical_space_->DoubleControlPointSlice(Dimension{dimension}, last);
    for (int index_point_slice_to_adjust = last; index_point_slice_to_adjust >= first; --index_point_slice_to_adjust) {
      util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> handler(this->GetNumberOfPointsPerDirection());
      --handler;
      handler.SetCurrentIndexForDimension(index_point_slice_to_adjust, Dimension{dimension});
      for (int i = handler.GetNumberOfMultiIndicesForCollapsedDimension(Dimension{dimension}) - 1; i >= 0;
           --i, handler.SubtractWithConstantDimension(Dimension{dimension})) {
        int current_upper_point_index = handler.GetCurrent1DIndex();
        int current_lower_point_index = handler.GetCurrent1DIndex() - handler.GetCurrentSliceSize(Dimension{dimension});
        double current_scaling = scaling[index_point_slice_to_adjust - first];
        auto new_control_point =
            ((1 - current_scaling) * physical_space_->GetHomogeneousControlPoint(current_lower_point_index) +
            current_scaling * physical_space_->GetHomogeneousControlPoint(current_upper_point_index));
        physical_space_->SetHomogeneousControlPoint(handler.GetCurrentIndex(), new_control_point);
      }
    }
  }

  bool RemoveControlPoints(std::vector<double> scaling, int first, int last, int dimension, double tolerance) override {
    int off = first - 1, i = first, j = last;
    std::vector<double> temp_w = GetTempNewWeights(scaling, off, last, i, j, dimension);
    std::vector<double> temp = GetTempNewControlPoints(scaling, temp_w, off, last, i, j, dimension);
    auto diff = static_cast<int>(ceil((j - i) / 2.0));
    i += diff, j -= diff;
    if (!IsKnotRemovable(scaling[i - off - 1], temp, temp_w, tolerance, i, j, off, dimension)) {
      return false;
    }
    SetNewControlPoints(temp, last, i - off, off, dimension);
    SetNewWeights(temp_w, last, i - off, off, dimension);
    physical_space_->RemoveControlPoints(
        this->GetTotalNumberOfControlPoints() / this->GetNumberOfPointsPerDirection()[dimension]);
    // TODO(all): Find a better solution than this one.
    physical_space_->DecrementNumberOfPoints(Dimension{dimension});
    return true;
  }

  std::array<std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>>, 2> SudivideSpline(ParametricCoordinate param_coord,
      int dimension) {
    this->InsertKnot(param_coord, dimension,
                     this->GetDegree(dimension).Get() + 1
                         - this->GetKnotVector(dimension)->GetMultiplicity(param_coord).Get());
    std::array<baf::KnotVectors<PARAMETRIC_DIMENSIONALITY>, 2>
        new_knot_vectors = this->parameter_space_->GetDividedKnotVectors(param_coord, dimension);
    std::array<Degree, PARAMETRIC_DIMENSIONALITY> degrees{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      degrees[i] = this->GetDegree(i);
    }
    std::array<std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>>, 2> subdivided_splines;
    int first = 0;
    for (int i = 0; i < 2; ++i) {
      int length = new_knot_vectors[i][dimension]->GetNumberOfKnots() - degrees[dimension].Get() - 1;
      std::vector<spl::ControlPoint> points = physical_space_->SplitControlPoints(Dimension{dimension}, first, length);
      std::vector<Weight> weights = physical_space_->SplitWeights(Dimension{dimension}, first, length);
      spl::NURBS<PARAMETRIC_DIMENSIONALITY> spline(new_knot_vectors[i], degrees, points, weights);
      subdivided_splines[i] = std::make_shared<spl::NURBS<PARAMETRIC_DIMENSIONALITY>>(spline);
      first = length;
    }
    return subdivided_splines;
  }

 private:
  std::shared_ptr<spl::PhysicalSpace<PARAMETRIC_DIMENSIONALITY>> GetPhysicalSpace() const override {
    return physical_space_;
  }

  double GetEvaluatedControlPoint(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                                  std::array<int, PARAMETRIC_DIMENSIONALITY> indices,
                                  int dimension) const override {
    if (this->GetPointDim() == dimension) {
      return this->parameter_space_->GetBasisFunctions(indices, param_coord)
          * physical_space_->GetHomogeneousControlPoint(indices)[Dimension{dimension}];
    }
    return this->parameter_space_->GetBasisFunctions(indices, param_coord)
        * physical_space_->GetHomogeneousControlPoint(indices)[Dimension{dimension}]
        / GetEvaluatedWeightSum(param_coord);
  }

  double GetEvaluatedDerivativeControlPoint(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                                            std::array<int, PARAMETRIC_DIMENSIONALITY> derivative,
                                            std::array<int, PARAMETRIC_DIMENSIONALITY> indices,
                                            int dimension) const override {
    return GetRationalBasisFunctionDerivative(param_coord, derivative, indices, dimension)
        * physical_space_->GetControlPoint(indices)[Dimension{dimension}];
  }

  double GetRationalBasisFunctionDerivative(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                                            std::array<int, PARAMETRIC_DIMENSIONALITY> derivative,
                                            std::array<int, PARAMETRIC_DIMENSIONALITY> indices,
                                            int dimension) const {
    if (derivative == std::array<int, PARAMETRIC_DIMENSIONALITY>{0}) {
      return this->parameter_space_->GetBasisFunctions(indices, param_coord)
          * physical_space_->GetWeight(indices).Get() / GetEvaluatedWeightSum(param_coord);
    }
    return (GetEvaluatedDerivativeWeight(param_coord, derivative, indices)
        - GetDerivativesSum(param_coord, derivative, indices, dimension))
        / GetEvaluatedWeightSum(param_coord);
  }

  double GetDerivativesSum(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                           std::array<int, PARAMETRIC_DIMENSIONALITY> derivative,
                           std::array<int, PARAMETRIC_DIMENSIONALITY> indices,
                           int dimension) const {
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> derivativeHandler(GetDerivativeHandler(derivative));
    derivativeHandler++;
    double sum = 0;
    for (int i = 1; i < derivativeHandler.GetNumberOfTotalMultiIndices(); i++, derivativeHandler++) {
      sum += binomialCoefficient(derivative, derivativeHandler.GetCurrentIndex())
          * GetRationalBasisFunctionDerivative(param_coord,
                                               derivativeHandler.GetComplementaryIndex(),
                                               indices,
                                               dimension)
          * GetEvaluatedDerivativeWeightSum(param_coord, derivativeHandler.GetCurrentIndex());
    }
    return sum;
  }

  double GetEvaluatedWeightSum(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord) const {
    auto first_non_zero = this->GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> basisFunctionHandler(
        this->GetNumberOfBasisFunctionsToEvaluate());
    double sum = 0;
    for (int i = 0; i < basisFunctionHandler.GetNumberOfTotalMultiIndices(); ++i, basisFunctionHandler++) {
      auto indices = basisFunctionHandler.GetCurrentIndex();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<>());
      sum += GetEvaluatedWeight(param_coord, indices);
    }
    return sum;
  }

  double GetEvaluatedDerivativeWeightSum(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                                         std::array<int, PARAMETRIC_DIMENSIONALITY> derivative) const {
    auto first_non_zero = this->GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> basisFunctionHandler(
        this->GetNumberOfBasisFunctionsToEvaluate());
    double sum = 0;
    for (int i = 0; i < basisFunctionHandler.GetNumberOfTotalMultiIndices(); ++i, basisFunctionHandler++) {
      auto indices = basisFunctionHandler.GetCurrentIndex();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<>());
      sum += GetEvaluatedDerivativeWeight(param_coord, derivative, indices);
    }
    return sum;
  }

  double GetEvaluatedWeight(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
      std::array<int, PARAMETRIC_DIMENSIONALITY> indices) const {
    return this->parameter_space_->GetBasisFunctions(indices, param_coord) * physical_space_->GetWeight(indices).Get();
  }

  double GetEvaluatedDerivativeWeight(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                                      std::array<int, PARAMETRIC_DIMENSIONALITY> derivative,
                                      std::array<int, PARAMETRIC_DIMENSIONALITY> indices) const {
    return this->parameter_space_->GetBasisFunctionDerivatives(indices, param_coord, derivative)
        * physical_space_->GetWeight(indices).Get();
  }

  util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> GetDerivativeHandler(
      const std::array<int, PARAMETRIC_DIMENSIONALITY> &derivative) const {
    std::array<int, PARAMETRIC_DIMENSIONALITY> derivative_length{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      derivative_length[i] = derivative[i] + 1;
    }
    return util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(derivative_length);
  }

  int binomialCoefficient(int number, int subset) const {
    if (subset == 0 || subset == number)
      return 1;
    return binomialCoefficient(number - 1, subset - 1) + binomialCoefficient(number - 1, subset);
  }

  int binomialCoefficient(std::array<int, PARAMETRIC_DIMENSIONALITY> numbers,
      std::array<int, PARAMETRIC_DIMENSIONALITY> subsets) const {
    int bc = 1;
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      bc *= binomialCoefficient(numbers[i], subsets[i]);
    }
    return bc;
  }

  void SetNewControlPoints(const std::vector<double> &temp, int last, int ii, int off, int dimension) {
    std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = this->GetNumberOfPointsPerDirection();
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
    for (int m = 0; m < point_handler.GetNumberOfTotalMultiIndices(); ++m, ++point_handler) {
      int k = point_handler[Dimension{dimension}];
      if (k - off >= 1 && k - off != ii && k < last + 2) {
        int index =
            (point_handler.CollapseDimension(Dimension{dimension}) * (last - off + 2) + k - off) * this->GetPointDim();
        std::vector<double> coordinates(temp.begin() + index, temp.begin() + index + this->GetPointDim());
        auto indices = point_handler.GetCurrentIndex();
        indices[dimension] = k - off < ii ? k : k - 1;
        physical_space_->SetControlPoint(indices, spl::ControlPoint(coordinates), Dimension{dimension},
                                         util::numeric_operations::decrement<int>);
      }
      if ((k <= off && k - off < 1) || (k >= last + 1 && k < this->GetNumberOfPointsPerDirection()[dimension])) {
        auto indices = point_handler.GetCurrentIndex();
        indices[dimension] = k <= off ? k : k - 1;
        physical_space_->SetControlPoint(indices, this->GetControlPoint(point_handler.GetCurrentIndex()),
                                         Dimension{dimension}, util::numeric_operations::decrement<int>);
      }
    }
  }

  void SetNewWeights(const std::vector<double> &temp, int last, int ii, int off, int dimension) {
    std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = this->GetNumberOfPointsPerDirection();
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
    for (int m = 0; m < point_handler.GetNumberOfTotalMultiIndices(); ++m, ++point_handler) {
      int k = point_handler[Dimension{dimension}];
      if (k - off >= 1 && k - off != ii && k < last + 2) {
        int index = point_handler.CollapseDimension(Dimension{dimension}) * (last - off + 2) + k - off;
        auto indices = point_handler.GetCurrentIndex();
        indices[dimension] = k - off < ii ? k : k - 1;
        physical_space_->SetWeight(indices, Weight{temp[index]}, Dimension{dimension},
                                   util::numeric_operations::decrement<int>);
      }
      if ((k <= off && k - off < 1) || (k >= last + 1 && k < this->GetNumberOfPointsPerDirection()[dimension])) {
        auto indices = point_handler.GetCurrentIndex();
        indices[dimension] = k <= off ? k : k - 1;
        physical_space_->SetWeight(indices, Weight{this->GetWeight(point_handler.GetCurrentIndex())},
                                   Dimension{dimension}, util::numeric_operations::decrement<int>);
      }
    }
  }

  std::vector<double> GetTempNewControlPoints(const std::vector<double> &scaling, const std::vector<double> &temp_w,
                                              int off, int last, int i, int j, int dimension) const {
    std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = this->GetNumberOfPointsPerDirection();
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
    int new_points = this->GetTotalNumberOfControlPoints() / this->GetNumberOfPointsPerDirection()[dimension];
    std::vector<double> temp(new_points * this->GetPointDim() * (last - off + 2), 0);
    std::shared_ptr<std::vector<double>> temp_ptr = std::make_shared<std::vector<double>>(temp);
    for (int l = 0; l < point_handler.GetNumberOfTotalMultiIndices(); ++l, ++point_handler) {
      if (point_handler[Dimension{dimension}] == off || point_handler[Dimension{dimension}] == last + 1) {
        int diff = point_handler[Dimension{dimension}] == off ? off : last + 1;
        SetTempNewControlPoint(point_handler, temp_ptr, temp_w, 1, diff, off, last, dimension, 0);
      }
    }
    for (; j - i > 0; ++i, --j) {
      point_handler.SetCurrentIndex({0});
      for (int l = 0; l < point_handler.GetNumberOfTotalMultiIndices(); ++l, ++point_handler) {
        if (point_handler[Dimension{dimension}] == i) {
          SetTempNewControlPoint(point_handler, temp_ptr, temp_w, scaling[i - off - 1], i, off, last, dimension, -1);
        }
        if (point_handler[Dimension{dimension}] == j) {
          SetTempNewControlPoint(point_handler, temp_ptr, temp_w, 1 - scaling[j - off - 1], j, off, last, dimension, 1);
        }
      }
    }
    return *temp_ptr;
  }

  void SetTempNewControlPoint(const util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> &point_handler,
                              const std::shared_ptr<std::vector<double>> &temp_ptr, const std::vector<double> &temp_w,
                              double alpha, int x, int off, int last, int dimension, int shift) const {
    int index = point_handler.CollapseDimension(Dimension{dimension}) * (last - off + 2);
    for (int k = 0; k < this->GetPointDim(); ++k) {
      (*temp_ptr)[(index + x - off) * this->GetPointDim() + k] =
          (GetHomogeneousControlPoint(point_handler.GetCurrentIndex(), k)
              - (1 - alpha) * (*temp_ptr)[(index + x - off + shift) * this->GetPointDim() + k]
                  * temp_w[index + x - off + shift]) / alpha / temp_w[index + x - off];
    }
  }

  std::vector<double> GetTempNewWeights(const std::vector<double> &scaling, int off, int last,
                                        int i, int j, int dimension) const {
    std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = this->GetNumberOfPointsPerDirection();
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
    int new_control_points = this->GetTotalNumberOfControlPoints() / this->GetNumberOfPointsPerDirection()[dimension];
    std::vector<double> temp_w(new_control_points * (last - off + 2), 0);
    std::shared_ptr<std::vector<double>> temp_w_ptr = std::make_shared<std::vector<double>>(temp_w);
    for (int l = 0; l < point_handler.GetNumberOfTotalMultiIndices(); ++l, ++point_handler) {
      if (point_handler[Dimension{dimension}] == off || point_handler[Dimension{dimension}] == last + 1) {
        int diff = point_handler[Dimension{dimension}] == off ? off : last + 1;
        SetTempNewWeight(point_handler, temp_w_ptr, 1, diff, off, last, dimension, 0);
      }
    }
    for (; j - i > 0; ++i, --j) {
      point_handler.SetCurrentIndex({0});
      for (int l = 0; l < point_handler.GetNumberOfTotalMultiIndices(); ++l, ++point_handler) {
        if (point_handler[Dimension{dimension}] == i) {
          SetTempNewWeight(point_handler, temp_w_ptr, scaling[i - off - 1], i, off, last, dimension, -1);
        }
        if (point_handler[Dimension{dimension}] == j) {
          SetTempNewWeight(point_handler, temp_w_ptr, 1 - scaling[j - off - 1], j, off, last, dimension, 1);
        }
      }
    }
    return *temp_w_ptr;
  }

  void SetTempNewWeight(const util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> &point_handler,
                        const std::shared_ptr<std::vector<double>> &temp_w_ptr,
                        double alpha, int x, int off, int last, int dimension, int shift) const {
    int index = point_handler.CollapseDimension(Dimension{dimension}) * (last - off + 2);
    (*temp_w_ptr)[index + x - off] = (this->GetWeight(point_handler.GetCurrentIndex())
        - (1 - alpha) * (*temp_w_ptr)[index + x - off + shift]) / alpha;
  }

  bool IsKnotRemovable(double alfi, const std::vector<double> &temp, const std::vector<double> &temp_w,
                       double tolerance, int i, int j, int off, int dimension) const {
    auto maxdist = physical_space_->GetMaximumDistanceFromOrigin();
    auto minw = physical_space_->GetMinimumWeight();
    tolerance = tolerance * minw / (1 + maxdist);

    std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = this->GetNumberOfPointsPerDirection();
    point_handler_length[dimension] = 0;
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
    int new_points = this->GetTotalNumberOfControlPoints() / this->GetNumberOfPointsPerDirection()[dimension];
    size_t temp_length = temp.size() / new_points;
    size_t temp_w_length = temp_w.size() / new_points;
    for (int l = 0; l < new_points; ++l, ++point_handler) {
      size_t offset = l * temp_length;
      size_t w_offset = l * temp_w_length;
      std::vector<double> temp1(this->GetPointDim() + 1, temp_w[w_offset + i - off - 1]);
      std::vector<double> temp2(this->GetPointDim() + 1, temp_w[w_offset + j - off + 1]);
      for (int k = 0; k < this->GetPointDim(); ++k) {
        temp1[k] = temp[offset + (i - off - 1) * this->GetPointDim() + k] * temp_w[w_offset + i - off - 1];
        temp2[k] = temp[offset + (j - off + 1) * this->GetPointDim() + k] * temp_w[w_offset + j - off + 1];
      }
      if (util::vector_utils::ComputeDistance(temp1, temp2) > tolerance) {
        auto indices = point_handler.GetCurrentIndex();
        indices[dimension] = i;
        for (int k = 0; k < this->GetPointDim(); ++k) {
          temp1[k] = physical_space_->GetHomogeneousControlPoint(indices)[Dimension{k}];
          temp2[k] = alfi * temp[offset + (i - off + 1) * this->GetPointDim() + k] * temp_w[w_offset + i - off + 1]
              + (1 - alfi) * temp[offset + (i - off - 1) * this->GetPointDim() + k] * temp_w[w_offset + i - off - 1];
        }
        temp1[this->GetPointDim()] = physical_space_->GetWeight(indices).Get();
        temp2[this->GetPointDim()] =
            alfi * temp_w[w_offset + i - off + 1] + (1 - alfi) * temp_w[w_offset + i - off - 1];
        if (util::vector_utils::ComputeDistance(temp1, temp2) > tolerance) {
          return false;
        }
      }
    }
    return true;
  }

  void SetNewControlPoint(spl::ControlPoint control_point, double weight,
      std::array<int, PARAMETRIC_DIMENSIONALITY> indices) override {
    physical_space_->SetWeightedControlPoint(indices, control_point, Weight{weight});
  }

  std::shared_ptr<WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>> physical_space_;
};
}  // namespace splinelib::src::spl

#endif  // SRC_SPL_NURBS_H_
