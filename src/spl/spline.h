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
#include <iostream>
#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <sstream>
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
    old_ = true;
  }
  explicit Spline(std::shared_ptr<ParameterSpace < DIM>>
parameter_space) {
    parameter_space_ = parameter_space;
    old_ = true;
  }
  Spline(const Spline<DIM> &spline) {
    ParameterSpace<DIM> parameter_space(*spline.parameter_space_);
    parameter_space_ = std::make_shared<ParameterSpace<DIM>>(parameter_space);
    old_ = true;
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

  void ElevateDegree(int dimension) {
    old_ = 0;
    std::vector<double> alpha = ComputeBezierDegreeElevationCoeffients(dimension);
    auto diff = ProduceBezierSegments(dimension);
    std::vector<std::vector<baf::ControlPoint>> new_bez_cps(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 1);
    for (int i = static_cast<int>(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2); i >= 0; --i) {
      new_bez_cps[i] = DegreeElevateBezierSegment(GetBezierSegment(dimension, i), alpha, dimension);
    }
    InsertKnot(GetKnotVector(dimension)->GetLastKnot(), dimension);
    parameter_space_->InsertKnot(GetKnotVector(dimension)->GetKnot(0), dimension);
    for (int i = static_cast<int>(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2); i >= 0; --i) {
      SetNewBezierSegmentControlPoints(new_bez_cps[i], dimension, i);
    }
    parameter_space_->ElevateDegree(dimension);
//    RemoveBezierKnots(diff, dimension);
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
    std::vector<int> diff(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2, GetDegree(dimension).get());
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
      RemoveKnot(GetKnotVector(dimension)->GetKnot(position),
                 dimension, 100 /*util::NumericSettings<double>::kEpsilon()*/, static_cast<size_t>(diff[i - 1] - 1));
    }
  }

  std::vector<double> GetBezierSegment(int dimension, int segment) const {
    int first = segment * (GetDegree(dimension).get() + 1);
    auto points_per_direction = GetPointsPerDirection();
    util::MultiIndexHandler<DIM> point_handler(points_per_direction);
    int width = GetDegree(dimension).get() + 1;
    int length = GetNumberOfControlPoints() / points_per_direction[dimension];
    int segment_width = GetPointDim() * (GetDegree(dimension).get() + 1);
    auto p = GetPointDim() * (GetDegree(dimension).get() + 1) * length;
    std::vector<double> bezier_cps(GetPointDim() * (GetDegree(dimension).get() + 1) * length);
    if (!old_) {
      for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
        auto a = point_handler.Get1DIndex();
        auto b1 = point_handler[0];
        auto b2 = point_handler[1];
        auto c = point_handler[dimension]; //point_handler.Get1DIndex() % points_per_direction[dimension];
        auto
            d = c >= segment * (GetDegree(dimension).get() + 1) && c < (segment + 1) * (GetDegree(dimension).get() + 1);
        auto e = point_handler.Get1DIndex() / points_per_direction[dimension];
        if (d) {
          for (int j = 0; j < GetPointDim(); ++j) {
            auto f = e * segment_width + (c % (GetDegree(dimension).get() + 1)) * GetPointDim() + j;
//            auto p = (c % width) * length + (point_handler.Get1DIndex() % c);
//            auto p1 = c % width;
//            auto p2 = point_handler.Get1DIndex() % c;
//            auto p3 = p1 * length;
//            auto p4 = p2 + p3;
            auto r = point_handler.ExtractDimension(dimension);
            auto q = c % width + point_handler.ExtractDimension(dimension) * width;
            auto qpd = q * GetPointDim() + j;
            auto g = GetControlPoint(point_handler.GetIndices(), j);
//            bezier_cps[e * segment_width + c * GetPointDim() + j] = GetControlPoint(point_handler.GetIndices(), j);
//            bezier_cps.emplace_back(GetControlPoint(point_handler.GetIndices(), j));
            bezier_cps[qpd] = GetControlPoint(point_handler.GetIndices(), j);
          }
        }
      }
    } else {
//      std::vector<double> bezier_cps(static_cast<size_t>(GetPointDim() * (GetDegree(dimension).get() + 1)));
      for (int i = 0; i < GetDegree(dimension).get() + 1; ++i) {
        for (int j = 0; j < GetPointDim(); ++j) {
//          bezier_cps[i * GetPointDim() + j] = GetControlPoint({first + i}, j);
          bezier_cps.emplace_back(GetControlPoint({first + i}, j));
        }
      }
    }

    std::cout << std::endl << "Bezier control points:" << std::endl;
    std::cout << bezier_cps.size() << std::endl;
    for (const auto p : bezier_cps) {
      std::cout << p << "  ";
    }
    std::cout << std::endl;

    return bezier_cps;
  }

  std::vector<baf::ControlPoint> DegreeElevateBezierSegment(std::vector<double> bezier_cps,
                                                            std::vector<double> alpha,
                                                            int dimension) const {
    auto points_per_direction = GetPointsPerDirection();
    int width = GetDegree(dimension).get() + 1;
    int length = GetNumberOfControlPoints() / points_per_direction[dimension];
    int segment_width = GetPointDim() * (GetDegree(dimension).get() + 1);
    util::MultiIndexHandler<2> point_handler({width, length});

    std::vector<baf::ControlPoint> cps;
    std::vector<double> coord(static_cast<size_t>(GetPointDim()));
    if (!old_) {

      for (int k = 0; k < point_handler.Get1DLength(); ++k, ++point_handler) {
        auto i = point_handler[0];
        auto p = point_handler.Get1DIndex(); // point_handler.ExtractDimension(dimension);
        for (int j = 0; j < GetPointDim(); ++j) {
          auto i1 = p * GetPointDim() + j;
          auto i2 = (p - 1) * GetPointDim() + j;
          auto c1 = bezier_cps[i1];
          auto c2 = bezier_cps[i2];
          auto e1 = (1 - alpha[i]);
          auto e2 = alpha[i];
          auto n = e1 * c1 + e2 * c2;
          coord[j] =
              (1 - alpha[i]) * bezier_cps[p * GetPointDim() + j] + alpha[i] * bezier_cps[(p - 1) * GetPointDim() + j];
        }
//        if (i >= 0 && i < GetDegree(dimension).get() + 1) {
//          for (int j = 0; j < GetPointDim(); ++j) {
//            coord[j] =
//                (1 - alpha[i]) * bezier_cps[p * GetPointDim() + j] + alpha[i] * bezier_cps[(p - 1) * GetPointDim() + j];
//          }
        cps.emplace_back(baf::ControlPoint(coord));
//        }
        if (point_handler[0] == GetDegree(dimension).get()) {
          for (int j = 0; j < GetPointDim(); ++j) {
            auto a = point_handler[0] * segment_width + GetDegree(dimension).get() * GetPointDim() + j;
            auto b = point_handler.Get1DIndex();
            coord[j] = bezier_cps[b * GetPointDim() + j];
          }
          cps.emplace_back(baf::ControlPoint(coord));
        }
      }
    } else {
      for (int i = 0; i < GetDegree(dimension).get() + 1; ++i) {
        for (int j = 0; j < GetPointDim(); ++j) {
          coord[j] =
              (1 - alpha[i]) * bezier_cps[i * GetPointDim() + j] + alpha[i] * bezier_cps[(i - 1) * GetPointDim() + j];
        }
        cps.emplace_back(baf::ControlPoint(coord));
      }
      for (int j = 0; j < GetPointDim(); ++j) {
        coord[j] = bezier_cps[GetDegree(dimension).get() * GetPointDim() + j];
      }
      cps.emplace_back(baf::ControlPoint(coord));
    }

    std::cout << std::endl << "degree elevated bezier segment: " << std::endl;
    std::cout << cps.size() << "  " << cps.size() * GetPointDim() << std::endl;
    for (const auto &y : cps) {
      for (int i = 0; i < GetPointDim(); ++i) {
        std::cout << y.GetValue(i) << "  ";
      }
    }
    std::cout << std::endl;

    return cps;
  }

  void SetNewBezierSegmentControlPoints(std::vector<baf::ControlPoint> new_bezier_cps, int dimension, int segment) {
    int first = segment * (GetDegree(dimension).get() + 1);
    auto points_per_direction = GetPointsPerDirection();
    util::MultiIndexHandler<DIM> point_handler(points_per_direction);
    int width = GetDegree(dimension).get() + 1;
    int length = GetNumberOfControlPoints() / points_per_direction[dimension];
    int segment_width = GetPointDim() * (GetDegree(dimension).get() + 1);
    for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
      auto c = point_handler[dimension]; // point_handler.Get1DIndex() % points_per_direction[dimension];
      auto
          d = c > segment * (GetDegree(dimension).get() + 1) && c <= (segment + 1) * (GetDegree(dimension).get() + 1);
      auto e = point_handler.Get1DIndex() / points_per_direction[dimension];
      auto f = new_bezier_cps.size();

      auto g = point_handler[0];
      auto h = point_handler[1];
      if (d) {
        auto m = e * width + c;
        auto n = e * (width + 1) + (c - 1) % (GetDegree(dimension).get() + 1) + 1;
        auto p1 = (c - 1) % width;
        auto p2 = point_handler.ExtractDimension(dimension);
        auto p3 = p2 * (width + 1);
        auto p4 = p1 + p3 + 1;
        auto q1 = new_bezier_cps[p4].GetValue(0);
        auto q2 = new_bezier_cps[p4].GetValue(1);
        auto s = new_bezier_cps[n].GetValue(0);
        auto t = new_bezier_cps[n].GetValue(1);
//        GetPhysicalSpace()->SetControlPoint(point_handler.GetIndices(), new_bezier_cps[n]);
        GetPhysicalSpace()->SetControlPoint(point_handler.GetIndices(), new_bezier_cps[p4]);
      }
    }


//    for (int i = 0; i < GetDegree(dimension).get() + 1; ++i) {
//      GetPhysicalSpace()->SetControlPoint({i + first}, new_bezier_cps[i]);
//    }
  }

  std::shared_ptr<ParameterSpace < DIM>> parameter_space_;
  bool old_;
};
}  //  namespace spl

#endif  //  SRC_SPL_SPLINE_H_
