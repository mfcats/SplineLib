/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_SPLINE_H_
#define SRC_SPL_SPLINE_H_

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <sstream>
#include <utility>
#include <vector>

#include "control_point.h"
#include "knot_vector.h"
#include "multi_index_handler.h"
#include "parameter_space.h"
#include "physical_space.h"

namespace spl {
template<int DIM>
class Spline {
 public:
  virtual ~Spline() = default;
  Spline() = default;
  Spline(KnotVectors<DIM> knot_vector, std::array<Degree, DIM> degree) {
    parameter_space_ = std::make_shared<ParameterSpace<DIM>>(ParameterSpace<DIM>(knot_vector, degree));
  }
  explicit Spline(std::shared_ptr<ParameterSpace<DIM>> parameter_space) {
    parameter_space_ = parameter_space;
  }
  Spline(const Spline<DIM> &spline) {
    ParameterSpace<DIM> parameter_space(*spline.parameter_space_);
    parameter_space_ = std::make_shared<ParameterSpace<DIM>>(parameter_space);
  }

  virtual std::vector<double> Evaluate(std::array<ParamCoord, DIM> param_coord,
                                       const std::vector<int> &dimensions) const {
    this->ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);

    auto first_non_zero = GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<DIM> basisFunctionHandler(this->GetNumberOfBasisFunctionsToEvaluate());
    std::vector<double> evaluated_point(dimensions.size(), 0);

    for (int i = 0; i < basisFunctionHandler.Get1DLength(); ++i, basisFunctionHandler++) {
      auto indices = basisFunctionHandler.GetIndices();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<double>());
      for (auto j = 0u; j < dimensions.size(); ++j) {
        evaluated_point[j] += GetEvaluatedControlPoint(param_coord, indices, dimensions[j]);
      }
    }
    return evaluated_point;
  }

  virtual std::vector<double> EvaluateDerivative(std::array<ParamCoord, DIM> param_coord,
                                                 const std::vector<int> &dimensions,
                                                 std::array<int, DIM> derivative) const {
    this->ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);

    auto first_non_zero = this->GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<DIM> basisFunctionHandler(this->GetNumberOfBasisFunctionsToEvaluate());
    std::vector<double> evaluated_point(dimensions.size(), 0);

    for (int i = 0; i < basisFunctionHandler.Get1DLength(); ++i, basisFunctionHandler++) {
      auto indices = basisFunctionHandler.GetIndices();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<double>());
      for (auto j = 0u; j < dimensions.size(); ++j) {
        evaluated_point[j] += GetEvaluatedDerivativeControlPoint(param_coord, derivative, indices, dimensions[j]);
      }
    }
    return evaluated_point;
  }

  std::vector<double> EvaluateAllNonZeroBasisFunctions(int direction, ParamCoord param_coord) const {
    return parameter_space_->EvaluateAllNonZeroBasisFunctions(direction, param_coord);
  }

  std::vector<double> EvaluateAllNonZeroBasisFunctionDerivatives(int direction,
                                                                 ParamCoord param_coord,
                                                                 int derivative) const {
    return parameter_space_->EvaluateAllNonZeroBasisFunctionDerivatives(direction, param_coord, derivative);
  }

  bool AreGeometricallyEqual(const spl::Spline<DIM> &rhs,
                             double tolerance = util::NumericSettings<double>::kEpsilon()) const {
    double number = ceil(pow(100, 1.0 / DIM));
    std::array<int, DIM> points = std::array<int, DIM>();
    std::vector<int> dimensions;
    for (int i = 0; i < DIM; ++i) {
      points[i] = number + 1;
    }
    for (int i = 0; i < GetPointDim(); ++i) {
      dimensions.push_back(i);
    }
    util::MultiIndexHandler<DIM> point_handler(points);
    for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
      std::array<ParamCoord, DIM> param_coord;
      for (int dim = 0; dim < DIM; ++dim) {
        double span = GetKnotVector(dim)->GetLastKnot().get() - GetKnotVector(dim)->GetKnot(0).get();
        param_coord[dim] = ParamCoord(span / number * point_handler[dim]);
      }
      std::vector<double> evaluate_this = Evaluate(param_coord, dimensions);
      std::vector<double> evaluate_rhs = rhs.Evaluate(param_coord, dimensions);
      if (util::VectorUtils<double>::ComputeDistance(evaluate_this, evaluate_rhs) > tolerance) {
        return false;
      }
    }
    return true;
  }

  Degree GetDegree(int i) const {
    return parameter_space_->GetDegree(i);
  }

  std::shared_ptr<baf::KnotVector> GetKnotVector(int i) const {
    return parameter_space_->GetKnotVector(i);
  }

  double GetKnotVectorRange(int direction) const {
    return parameter_space_->GetKnotVectorRange(direction);
  }

  int GetNumberOfControlPoints() const {
    return GetPhysicalSpace()->GetNumberOfControlPoints();
  }

  std::array<int, DIM> GetPointsPerDirection() const {
    return GetPhysicalSpace()->GetPointsPerDirection();
  }

  int GetPointDim() const {
    return GetPhysicalSpace()->GetDimension();
  }

  double GetControlPoint(std::array<int, DIM> indices, int dimension) const {
    return GetPhysicalSpace()->GetControlPoint(indices).GetValue(dimension);
  }

  baf::ControlPoint GetControlPoint(std::array<int, DIM> indices) const {
    return GetPhysicalSpace()->GetControlPoint(indices);
  }

  double GetExpansion() const {
    return GetPhysicalSpace()->GetExpansion();
  }

  double GetWeight(std::array<int, DIM> indices) const {
    return GetPhysicalSpace()->GetWeight(indices);
  }

  void InsertKnot(ParamCoord knot, int dimension, size_t multiplicity = 1) {
    KnotSpan knot_span = GetKnotVector(dimension)->GetKnotSpan(knot);
    Degree degree = GetDegree(dimension);
    for (size_t i = 1; i <= multiplicity; ++i) {
      if (GetKnotVector(dimension)->IsLastKnot(knot)) {
        knot_span = knot_span + KnotSpan{1};
      }
      auto last = knot_span.get() - GetKnotVector(dimension)->GetMultiplicity(knot);
      auto first = knot_span.get() - degree.get() + i;
      std::vector<double> scaling;
      for (auto j = static_cast<size_t>(first); j <= last; ++j) {
        ParamCoord low_knot = GetKnotVector(dimension)->GetKnot(j);
        ParamCoord upper_knot = GetKnotVector(dimension)->GetKnot(j + degree.get() - i + 1);
        scaling.emplace_back((knot.get() - low_knot.get()) / (upper_knot.get() - low_knot.get()));
      }
      this->AdjustControlPoints(scaling, static_cast<int>(first), static_cast<int>(last), dimension);
    }
    for (size_t i = 0; i < multiplicity; ++i) {
      parameter_space_->InsertKnot(knot, dimension);
    }
  }

  void RefineKnots(std::vector<ParamCoord> new_knots, int dimension) {
    for (const auto &knot : new_knots) {
      this->InsertKnot(knot, dimension);
    }
  }

  size_t RemoveKnot(ParamCoord knot, int dimension, double tolerance, size_t multiplicity = 1) {
    size_t count = 0;
    for (; count < multiplicity; ++count) {
      KnotSpan knot_span = GetKnotVector(dimension)->GetKnotSpan(knot);
      Degree degree = GetDegree(dimension);
      auto last = knot_span.get() - GetKnotVector(dimension)->GetMultiplicity(knot);
      auto first = knot_span.get() - degree.get();
      std::vector<double> scaling;
      for (auto i = static_cast<size_t>(first); i <= last; ++i) {
        ParamCoord low_knot = GetKnotVector(dimension)->GetKnot(i);
        ParamCoord upper_knot = GetKnotVector(dimension)->GetKnot(i + degree.get() + 1);
        scaling.emplace_back((knot.get() - low_knot.get()) / (upper_knot.get() - low_knot.get()));
      }
      bool is_removed = this->RemoveControlPoints(scaling, first, static_cast<int>(last), dimension, tolerance);
      if (!is_removed) break;
      parameter_space_->RemoveKnot(knot, dimension);
    }
    return count;
  }

  virtual void ElevateDegree(int dimension) {
    std::vector<double> alpha = ComputeBezierDegreeElevationCoeffients(dimension);
    auto diff = ProduceBezierSegments(dimension);
    std::vector<std::vector<std::pair<baf::ControlPoint, double>>>
        new_bez_points(this->GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 1);
    for (int i = static_cast<int>(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2); i >= 0; --i) {
      new_bez_points[i] = DegreeElevateBezierSegment(GetBezierSegment(dimension, i), alpha, dimension);
    }

    InsertKnot(GetKnotVector(dimension)->GetLastKnot(), dimension);
    // parameter_space_->InsertKnot(GetKnotVector(dimension)->GetLastKnot(), dimension);
    parameter_space_->InsertKnot(GetKnotVector(dimension)->GetKnot(0), dimension);

    // TODO: CODE CHANGED HERE!
    size_t position = GetKnotVector(dimension)->GetNumberOfKnots() - GetDegree(dimension).get() - 3;
    for (size_t i = GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2; i > 0; --i, --position) {
      while (GetKnotVector(dimension)->GetKnot(position) == GetKnotVector(dimension)->GetKnot(position - 1)) {
        --position;
      }
      InsertKnot(GetKnotVector(dimension)->GetKnot(position), dimension, 1);
    }

    for (int i = static_cast<int>(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2); i >= 0; --i) {
      SetNewBezierSegmentControlPoints(new_bez_points[i], dimension, i);
    }
    parameter_space_->ElevateDegree(dimension);

    // TODO: CODE CHANGED HERE!
    std::transform(diff.begin(), diff.end(), diff.begin(), std::bind(std::plus<int>(), std::placeholders::_1, 1));

    RemoveBezierKnots(diff, dimension);
  }

  void ReduceDegree(int dimension) {

    // TODO: Try to move the computation of the degree reduction coefficients to a separate function!
    // Compute p degree elevation coefficients that would be used for the degree elevation from p-1 to p.
    // std::vector<double> alpha = ComputeBezierDegreeReductionCoefficients(dimension);

    // Number of Bézier segments obtained by the decomposition below.
    int num_bezier_segments = GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 1;

    // Insert knots to subdivide the spline into Bézier segments. The diff vector holds the number of times every
    // inner knot was inserted in order to reach a multiplicity of p.
    std::vector<int> diff = ProduceBezierSegments(dimension);

    // Data structure to hold the new control points for all Bézier segments. First index is the number of the Bézier
    // segment, second index is the number of the control point of the specific Bézier segment. The double associated
    // with the control point (using std::pair) is the new weight.
    std::vector<std::vector<std::pair<baf::ControlPoint, double>>> all_new_bezier_cps;

    // Loop through all Bézier segments.
    for (int i = 0; i < num_bezier_segments; ++i) {
      all_new_bezier_cps.emplace_back(DegreeReduceBezierSegment(GetBezierSegment(dimension, i), dimension));
    }

    // Assemble the computed control points into the physical space of the spline.
    std::array<int, 1> indices = {1};
    int num_cps_per_bezier_segment = GetDegree(dimension).get();
    for (int i = 0; i < num_bezier_segments; ++i) {

      // TODO: Modify function below to be able to use it for degree elevation and degree reduction.
      // SetNewBezierSegmentControlPoints(all_new_bezier_cps[i], dimension, i);
      // Loop through all new control points of the current Bézier segment.
      for (int j = 1; j < num_cps_per_bezier_segment; ++j) {
        SetNewBezierPoint(all_new_bezier_cps[i][j], indices);
        ++indices[0];
      }
    }

    // TODO: On every Bézier segment, the number of control points is decreased by one.
    for (int i = 0; i < num_bezier_segments; ++i) {
      this->GetPhysicalSpace()->RemoveControlPoints(1);
      this->GetPhysicalSpace()->DecrementNumberOfPoints(0);
    }

    // TODO: Account for multiplicity of p+1 at inner knot values (see ProduceBezierSegments function).
    parameter_space_->RemoveKnot(ParamCoord{0.3}, 0);
    parameter_space_->RemoveKnot(ParamCoord{0.6}, 0);

    // Remove one entry of the first and last knot.
    parameter_space_->RemoveKnot(GetKnotVector(dimension)->GetKnot(0), dimension);
    parameter_space_->RemoveKnot(GetKnotVector(dimension)->GetLastKnot(), dimension);

    // Reduce degree of dimension by one and create new basis functions.
    parameter_space_->ReduceDegree(dimension);

    // Remove the knots to get rid of the decomposition into bezier segments (remove one knot more than inserted).
    std::transform(diff.begin(), diff.end(), diff.begin(), std::bind(std::plus<int>(), std::placeholders::_1, 1));
    RemoveBezierKnots(diff, dimension);
  }

  std::vector<std::pair<baf::ControlPoint, double>> DegreeReduceBezierSegment(std::vector<double> bezier_cps,
      int dimension) const {

    // Number of control points per degree reduced Bézier segment.
    uint64_t num_bezier_cps = GetDegree(dimension).get();

    // TODO: Use weights of one for now. Needs to be corrected to handle NURBS correctly!
    double new_weight = 1.0;

    // The first and last control point of the Bézier segment remain unchanged.
    baf::ControlPoint first(std::vector<double>(bezier_cps.begin(), bezier_cps.begin() + GetPointDim()));
    baf::ControlPoint last(std::vector<double>(bezier_cps.end() - GetPointDim(), bezier_cps.end()));

    // Vector of pairs to hold control point and corresponding weight.
    // TODO: Is there a cleaner way to allocate memory when creating the vector?
    std::vector<std::pair<baf::ControlPoint, double>> weighted_cps(num_bezier_cps, std::make_pair(first, new_weight));
    // weighted_cps[0] = std::make_pair(first, new_weight);
    weighted_cps[num_bezier_cps - 1] = std::make_pair(last, new_weight);

    // Vector to hold the physical coordinates of one control point.
    std::vector<double> coord(static_cast<size_t>(GetPointDim()));

    uint64_t r = (GetDegree(dimension).get() - 1) / 2;

    if (GetDegree(dimension).get() % 2 == 0) {
      // Loop through the new control points P_i for i = 1,...,r.
      for (uint64_t i = 1; i <= r; ++i) {

        // Compute the degree elevation coeffiecient alpha_i for the current control point.
        double alpha_i = double(i) / double(GetDegree(dimension).get());

        // Loop through the physical dimensions of the current new control point.
        // TODO: Control point times scalar should be implemented as an operator in the corresponding class!
        for (int j = 0; j < GetPointDim(); ++j) {

          // Compute the new control point P_i as a function of P_i-1, the old control point Q_i and the degree
          // elevation coefficient alpha_i.
          coord[j] =
              (bezier_cps[i * GetPointDim() + j] - alpha_i * weighted_cps[i - 1].first.GetValue(j)) / (1 - alpha_i);
        }
        weighted_cps[i] = std::make_pair(baf::ControlPoint(coord), new_weight);
      }

      // Loop though the new control points P_i for i = p-2,...,r+1.
      // Depending on the degree, p-2 can be greater or less than r+1. Therefore, the lower and upper bounds of the
      // loop are determined as follows. TODO: Correct?
      int lower = std::min(static_cast<int>(r + 1), GetDegree(dimension).get() - 2);
      int upper = std::max(static_cast<int>(r + 1), GetDegree(dimension).get() - 2);
      for (int i = upper; i >= lower; --i) {

        // Compute the degree elevation coefficient alpha_i for the current control point.
        double alpha_i = double(i + 1) / double(GetDegree(dimension).get());

        // Loop through the physical dimensions of the current new control point.
        for (int j = 0; j < GetPointDim(); ++j) {

          // Compute the new control point P_i as a function of P_i-1, the old control point Q_i and the degree
          // elevation coefficient alpha_i.
          coord[j] = (bezier_cps[(i + 1) * GetPointDim() + j]
              - (1 - alpha_i) * weighted_cps[i + 1].first.GetValue(j)) / alpha_i;
        }
        weighted_cps[i] = std::make_pair(baf::ControlPoint(coord), new_weight);
      }
    }

    if (GetDegree(dimension).get() % 2 != 0) {
      // Loop through the new control points P_i for i = 1,...,r-1.
      for (uint64_t i = 1; i <= r - 1; ++i) {

        // Compute the degree elevation coeffiecient alpha_i for the current control point.
        double alpha_i = double(i) / double(GetDegree(dimension).get());

        // Loop through the physical dimensions of the current new control point.
        // TODO: Control point times scalar should be implemented as an operator in the corresponding class!
        for (int j = 0; j < GetPointDim(); ++j) {

          // Compute the new control point P_i as a function of P_i-1, the old control point Q_i and the degree
          // elevation coefficient alpha_i.
          coord[j] =
              (bezier_cps[i * GetPointDim() + j] - alpha_i * weighted_cps[i - 1].first.GetValue(j)) / (1 - alpha_i);
        }
        weighted_cps[i] = std::make_pair(baf::ControlPoint(coord), new_weight);
      }

      // Loop though the new control points P_i for i = p-2,...,r+1.
      // Depending on the degree, p-2 can be greater or less than r+1. Therefore, the lower and upper bounds of the
      // loop are determined as follows. TODO: Correct?
      int lower = std::min(static_cast<int>(r + 1), GetDegree(dimension).get() - 2);
      int upper = std::max(static_cast<int>(r + 1), GetDegree(dimension).get() - 2);
      for (int i = upper; i >= lower; --i) {

        // Compute the degree elevation coefficient alpha_i for the current control point.
        double alpha_i = double(i + 1) / double(GetDegree(dimension).get());

        // Loop through the physical dimensions of the current new control point.
        for (int j = 0; j < GetPointDim(); ++j) {

          // Compute the new control point P_i as a function of P_i-1, the old control point Q_i and the degree
          // elevation coefficient alpha_i.
          coord[j] = (bezier_cps[(i + 1) * GetPointDim() + j]
              - (1 - alpha_i) * weighted_cps[i + 1].first.GetValue(j)) / alpha_i;
        }
        weighted_cps[i] = std::make_pair(baf::ControlPoint(coord), new_weight);
      }

      // Compute the control point in the center of the Bézier segment.
      double alpha_r_L = double(r) / double(GetDegree(dimension).get());
      double alpha_r_R = double(r + 1) / double(GetDegree(dimension).get());
      for (int j = 0; j < GetPointDim(); ++j) {
        double P_r_L = (bezier_cps[r * GetPointDim() + j] - alpha_r_L * weighted_cps[r - 1].first.GetValue(j))
            / (1 - alpha_r_L);
        double P_r_R = (bezier_cps[(r + 1) * GetPointDim() + j]
            - (1 - alpha_r_R) * weighted_cps[r + 1].first.GetValue(j)) / alpha_r_R;
        coord[j] = 0.5 * (P_r_L + P_r_R);
      }
      weighted_cps[r] = std::make_pair(baf::ControlPoint(coord), new_weight);
    }

    return weighted_cps;
  }

  virtual void AdjustControlPoints(std::vector<double> scaling, int first, int last, int dimension) = 0;
  virtual bool RemoveControlPoints(std::vector<double> scaling,
                                   int first, int last, int dimension, double tolerance) = 0;

 protected:
  void ThrowIfParametricCoordinateOutsideKnotVectorRange(std::array<ParamCoord, DIM> param_coord) const {
    parameter_space_->ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);
  }

  virtual double GetEvaluatedControlPoint(std::array<ParamCoord, DIM> param_coord,
                                          std::array<int, DIM> indices,
                                          int dimension) const = 0;

  virtual double GetEvaluatedDerivativeControlPoint(std::array<ParamCoord, DIM> param_coord,
                                                    std::array<int, DIM> derivative,
                                                    std::array<int, DIM> indices,
                                                    int dimension) const = 0;

  virtual std::shared_ptr<spl::PhysicalSpace<DIM>> GetPhysicalSpace() const = 0;

  std::array<int, DIM> GetArrayOfFirstNonZeroBasisFunctions(std::array<ParamCoord, DIM> param_coord) const {
    return parameter_space_->GetArrayOfFirstNonZeroBasisFunctions(param_coord);
  }

  std::array<int, DIM> GetNumberOfBasisFunctionsToEvaluate() const {
    std::array<int, DIM> total_length;
    for (int i = 0; i < DIM; ++i) {
      total_length[i] = parameter_space_->GetDegree(i).get() + 1;
    }
    return total_length;
  }

  std::vector<double> ComputeBezierDegreeElevationCoeffients(int dimension) const {
    std::vector<double> alpha;
    for (int i = 0; i < GetDegree(dimension).get() + 2; ++i) {
      alpha.push_back(static_cast<double>(i) / (GetDegree(dimension).get() + 1));
    }
    return alpha;
  }

  std::vector<int> ProduceBezierSegments(int dimension) {
    size_t position = GetKnotVector(dimension)->GetNumberOfKnots() - GetDegree(dimension).get() - 2;

    // TODO: In order to decompose the spline into Bézier segments, every inner knot has to be repeated only p times.
    // std::vector<int> diff(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2, GetDegree(dimension).get());
    std::vector<int> diff(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2, GetDegree(dimension).get() - 1);

    for (size_t i = GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2; i > 0; --i, --position) {
      while (GetKnotVector(dimension)->GetKnot(position) == GetKnotVector(dimension)->GetKnot(position - 1)) {
        --position;
        --diff[i - 1];
      }
      InsertKnot(GetKnotVector(dimension)->GetKnot(position), dimension, static_cast<size_t>(diff[i - 1]));
    }
    return diff;
  }

  void RemoveBezierKnots(std::vector<int> diff, int dimension) {
    size_t position = GetKnotVector(dimension)->GetNumberOfKnots() - GetDegree(dimension).get() - 2;
    for (size_t i = GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2; i > 0; --i, --position) {
      while (GetKnotVector(dimension)->GetKnot(position) == GetKnotVector(dimension)->GetKnot(position - 1)) {
        --position;
      }
      RemoveKnot(GetKnotVector(dimension)->GetKnot(position), dimension, util::NumericSettings<double>::kEpsilon(),
          static_cast<size_t>(diff[i - 1] - 1));
    }
  }

  virtual int GetBezierPointLength() const = 0;

  virtual std::vector<double> GetBezierSegment(int dimension, int segment) const {
    util::MultiIndexHandler<DIM> point_handler(GetPointsPerDirection());
    int width = GetDegree(dimension).get() + 1;
    int segment_length = GetNumberOfControlPoints() / GetPointsPerDirection()[dimension];
    int point_length = this->GetBezierPointLength();
    std::vector<double> bezier_cps(static_cast<size_t>(point_length * width * segment_length));
    for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
      auto index_in_dir = point_handler[dimension];

      // TODO: If the decomposition into Bézier segments is done correctly (insert knots up to multiplicity p), the line
      // below has to be modified to catch the right control points for all Bézier segments. As it is it will only
      // work for the first Bézier segment and will miss the first control point for the following Bézier segment.
      // if (index_in_dir >= segment * width && index_in_dir < (segment + 1) * width) {

      if (index_in_dir >= segment * width - segment && index_in_dir < (segment + 1) * width - segment) {
        double weight = this->GetWeight(point_handler.GetIndices());
        if (point_length == GetPointDim() + 1) {

          // TODO: See the comment above to understand the change here!
          // auto weight_index = (index_in_dir % width + point_handler.ExtractDimension(dimension) * width) * point_length
          //     + this->GetPointDim();
          auto weight_index = ((index_in_dir + segment) % width
              + point_handler.ExtractDimension(dimension) * width) * point_length + this->GetPointDim();

          bezier_cps[weight_index] = weight;
        }
        for (int j = 0; j < GetPointDim(); ++j) {

          // TODO: See the comment above to understand the change here!
          // auto index = (index_in_dir % width + point_handler.ExtractDimension(dimension) * width) * point_length + j;
          auto index = ((index_in_dir + segment) % width
              + point_handler.ExtractDimension(dimension) * width) * point_length + j;

          bezier_cps[index] = GetControlPoint(point_handler.GetIndices(), j) * weight;
        }
      }

    }
    return bezier_cps;
  }

  std::vector<std::pair<baf::ControlPoint, double>>
  DegreeElevateBezierSegment(std::vector<double> bezier_cps, std::vector<double> alpha, int dimension) const {
    int width = GetDegree(dimension).get() + 1;
    int segment_length = GetNumberOfControlPoints() / GetPointsPerDirection()[dimension];
    util::MultiIndexHandler<2> point_handler({width, segment_length});
    int point_length = this->GetBezierPointLength();
    std::vector<std::pair<baf::ControlPoint, double>> weighted_cps;
    std::vector<double> coord(static_cast<size_t>(GetPointDim()));
    for (int k = 0; k < point_handler.Get1DLength(); ++k, ++point_handler) {
      auto index = point_handler[0];
      auto pos = point_handler.Get1DIndex();
      double new_weight = point_length == GetPointDim() ? 1.0 :
                          (1 - alpha[index]) * bezier_cps[pos * point_length + this->GetPointDim()]
                              + alpha[index] * bezier_cps[(pos - 1) * point_length + this->GetPointDim()];
      for (int j = 0; j < GetPointDim(); ++j) {
        coord[j] = ((1 - alpha[index]) * bezier_cps[pos * point_length + j]
            + alpha[index] * bezier_cps[(pos - 1) * point_length + j]) / new_weight;
      }
      weighted_cps.emplace_back(std::make_pair(baf::ControlPoint(coord), new_weight));
      if (point_handler[0] == GetDegree(dimension).get()) {
        new_weight = point_length == GetPointDim() ? 1.0 : bezier_cps[pos * point_length + this->GetPointDim()];
        for (int j = 0; j < GetPointDim(); ++j) {
          coord[j] = bezier_cps[pos * point_length + j] / new_weight;
        }
        weighted_cps.emplace_back(std::make_pair(baf::ControlPoint(coord), new_weight));
      }
    }
    return weighted_cps;
  }

  virtual void SetNewBezierPoint(std::pair<baf::ControlPoint, double> new_bezier_point,
                                 std::array<int, DIM> indices) = 0;

  virtual void SetNewBezierSegmentControlPoints(
      const std::vector<std::pair<baf::ControlPoint, double>> &new_bezier_points, int dimension, int segment) {
    util::MultiIndexHandler<DIM> point_handler(GetPointsPerDirection());
    int width = GetDegree(dimension).get() + 1;
    for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
      auto index_in_dir = point_handler[dimension];
      if (index_in_dir > segment * width && index_in_dir <= (segment + 1) * width) {
        auto index = (index_in_dir - 1) % width + point_handler.ExtractDimension(dimension) * (width + 1) + 1;
        this->SetNewBezierPoint(new_bezier_points[index], point_handler.GetIndices());
      }
    }
  }

  std::shared_ptr<ParameterSpace<DIM>> parameter_space_;
};
}  // namespace spl

#endif  // SRC_SPL_SPLINE_H_
