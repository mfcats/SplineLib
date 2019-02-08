/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_PROJECTION_H_
#define SRC_SPL_PROJECTION_H_

#include <vector>

#include "b_spline.h"
#include "element_generator.h"
#include "knot_vector.h"
#include "vector_utils.h"

namespace spl {
template<int DIM>
class Projection {
 public:
  static std::array<ParamCoord, DIM> ProjectionOnSurface(const std::vector<double> &point_phys_coords,
                                                         std::shared_ptr<spl::Spline<DIM>> spline,
                                                         std::array<int, 2> scattering = {10, 10},
                                                         double tolerance = 0.00001) {
    int iteration = 0;
    bool converged = false;
    std::vector<int> dimensions = GetDimensions(point_phys_coords);
    std::array<ParamCoord, DIM> param_coords = FindInitialValue2D(scattering, point_phys_coords, spline, dimensions);
    while (!converged) {
      std::vector<double> derivative_dir1 = spline->EvaluateDerivative(param_coords, dimensions, {1, 0});
      std::vector<double> derivative_dir2 = spline->EvaluateDerivative(param_coords, dimensions, {0, 1});
      std::vector<double> current_point = spline->Evaluate(param_coords, dimensions);
      std::array<double, 2> delta = GetDelta2D(current_point, derivative_dir1, derivative_dir2, point_phys_coords);
      param_coords = UpdateProjectedPoint(param_coords, delta);
      param_coords = CheckParametricCoordinates(param_coords, spline);
      converged = UpdateConverged(delta, tolerance);
      if (++iteration > 1000) {
        break;
      }
    }
    return param_coords;
  }

  static std::array<ParamCoord, DIM> ProjectionOnCurve(const std::vector<double> &point_phys_coords,
                                                       std::shared_ptr<spl::Spline<DIM>> spline,
                                                       double tolerance = 0.0001) {
    int iteration = 0;
    bool converged = false;
    std::vector<int> dimensions = GetDimensions(point_phys_coords);
    std::array<ParamCoord, DIM> param_coords = FindInitialValue1D(point_phys_coords, spline, dimensions);
    while (!converged) {
      std::vector<double> derivative = spline->EvaluateDerivative(param_coords, dimensions, {1});
      std::vector<double> current_point = spline->Evaluate(param_coords, dimensions);
      std::array<double, DIM> delta = GetDelta1D(current_point, derivative, point_phys_coords);
      param_coords = UpdateProjectedPoint(param_coords, delta);
      param_coords = CheckParametricCoordinates(param_coords, spline);
      converged = UpdateConverged(delta, tolerance);
      if (++iteration > 1000) {
        break;
      }
    }
    return param_coords;
  }

 private:
  static std::vector<int> GetDimensions(const std::vector<double> &point_phys_coords) {
    std::vector<int> dimensions;
    for (auto i = 0u; i < point_phys_coords.size(); ++i) {
      dimensions.emplace_back(i);
    }
    return dimensions;
  }

  static std::array<ParamCoord, DIM> FindInitialValue1D(const std::vector<double> &point_phys_coords,
                                                        const std::shared_ptr<spl::Spline<DIM>> &spline,
                                                        const std::vector<int> &dimensions) {
    util::ElementGenerator<1> element_generator(spline);
    std::vector<util::Element> elements = element_generator.GetElementList(0);
    std::array<ParamCoord, DIM> paramCoords = {ParamCoord{0}};
    std::vector<double> splinePhysicalCoords = spline->Evaluate({ParamCoord{(0.5 * (
        elements[0].GetUpperBound() - elements[0].GetLowerBound()).get())}}, dimensions);
    double distance = util::VectorUtils<double>::ComputeDistance(point_phys_coords, splinePhysicalCoords);
    paramCoords[0] = ParamCoord{{0.5 * (elements[0].GetUpperBound() - elements[0].GetLowerBound()).get()}};
    for (auto i = 1u; i < elements.size(); ++i) {
      splinePhysicalCoords = spline->Evaluate({ParamCoord{0.5 * (
          elements[i].GetUpperBound() - elements[i].GetLowerBound()).get() + elements[i].GetLowerBound().get()}},
              dimensions);
      if (util::VectorUtils<double>::ComputeDistance(point_phys_coords, splinePhysicalCoords) < distance) {
        distance = util::VectorUtils<double>::ComputeDistance(point_phys_coords, splinePhysicalCoords);
        paramCoords[0] = ParamCoord{0.5 * (
                elements[i].GetUpperBound() - elements[i].GetLowerBound()).get() + elements[i].GetLowerBound().get()};
      }
    }
    return paramCoords;
  }

  static std::array<ParamCoord, 2> FindInitialValue2D(std::array<int, 2> scattering,
                                                      const std::vector<double> &point_phys_coords,
                                                      const std::shared_ptr<spl::Spline<2>> &spline,
                                                      const std::vector<int> &dimensions) {
    double first_knot1 = spline->GetKnotVector(0)->GetKnot(0).get();
    double last_knot1 = spline->GetKnotVector(0)->GetLastKnot().get();
    double first_knot2 = spline->GetKnotVector(1)->GetKnot(0).get();
    double last_knot2 = spline->GetKnotVector(1)->GetLastKnot().get();
    std::array<ParamCoord, 2> param_coords = {ParamCoord(first_knot1), ParamCoord(first_knot2)};
    std::vector<double> physical_coords = spline->Evaluate(param_coords, dimensions);
    double distance = util::VectorUtils<double>::ComputeDistance(point_phys_coords, physical_coords);
    for (int i = 1; i <= scattering[0]; ++i) {
      for (int j = 1; j <= scattering[1]; ++j) {
        ParamCoord coord1 = ParamCoord(i * (last_knot1 - first_knot1) / scattering[0]);
        ParamCoord coord2 = ParamCoord(j * (last_knot2 - first_knot2) / scattering[1]);
        physical_coords = spline->Evaluate({coord1, coord2}, dimensions);
        double current_dist = util::VectorUtils<double>::ComputeTwoNorm(util::VectorUtils<double>::ComputeDifference(
            point_phys_coords,
            physical_coords));
        if (current_dist < distance) {
          param_coords = {coord1, coord2};
          distance = current_dist;
        }
      }
    }
    return param_coords;
  }

  static std::array<ParamCoord, DIM> UpdateProjectedPoint(std::array<ParamCoord, DIM> param_coords,
                                                          const std::array<double, DIM> delta) {
    for (int i = 0; i < DIM; ++i) {
      param_coords[i] = param_coords[i] + ParamCoord{delta[i]};
    }
    return param_coords;
  }

  static std::array<ParamCoord, DIM> CheckParametricCoordinates(std::array<ParamCoord, DIM> param_coords,
                                                                const std::shared_ptr<spl::Spline<DIM>> &spline) {
    for (auto &coord : param_coords) {
      if (coord < spline->GetKnotVector(0)->GetKnot(0)) {
        coord = spline->GetKnotVector(0)->GetKnot(0);
      } else if (coord > spline->GetKnotVector(0)->GetLastKnot()) {
        coord = spline->GetKnotVector(0)->GetLastKnot();
      }
    }
    return param_coords;
  }

  static bool UpdateConverged(const std::array<double, DIM> &delta, double tolerance) {
    bool is_below_tolerance = true;
    for (const auto &i : delta) {
      if (std::abs(i) >= tolerance) {
        is_below_tolerance = false;
      }
    }
    return is_below_tolerance;
  }

  static std::vector<double> GetQ(const std::vector<double> &projection_point,
                                  const std::vector<double> &current_point,
                                  const std::vector<double> &derivative_dir1,
                                  const std::vector<double> &derivative_dir2) {
    std::vector<double> normal_vector = util::VectorUtils<double>::CrossProduct(derivative_dir1, derivative_dir2);
    double length = util::VectorUtils<double>::ComputeTwoNorm(normal_vector);
    std::vector<double> normed_normal_vector = util::VectorUtils<double>::ScaleVector(normal_vector, 1 / length);
    double distance = util::VectorUtils<double>::ComputeScalarProduct(current_point, normed_normal_vector)
        - util::VectorUtils<double>::ComputeScalarProduct(projection_point, normed_normal_vector);
    std::vector<double> scaled_normal_vector = util::VectorUtils<double>::ScaleVector(normed_normal_vector, -distance);
    std::vector<double> q = util::VectorUtils<double>::ComputeDifference(projection_point, scaled_normal_vector);
    return q;
  }

  static std::array<double, 1> GetDelta1D(const std::vector<double> &current_point,
                                          const std::vector<double> &derivative,
                                          const std::vector<double> &point_phys_coords) {
    std::vector<double> projectionVector =
        util::VectorUtils<double>::ComputeDifference(point_phys_coords, current_point);
    // This is only the first order algorithm. An implemented but non-working version of the second order algorithm
    // can be found in commit 2ed993e6dcef3d184b70640f6b9498efae52747a.
    return {util::VectorUtils<double>::ComputeScalarProduct(derivative, projectionVector)
                / util::VectorUtils<double>::ComputeScalarProduct(derivative, derivative)};
  }

  static std::array<double, 2> GetDelta2D(const std::vector<double> &current_point,
                                          const std::vector<double> &derivative_dir1,
                                          const std::vector<double> &derivative_dir2,
                                          const std::vector<double> &point_phys_coords) {
    std::vector<double> q = GetQ(point_phys_coords, current_point, derivative_dir1, derivative_dir2);
    std::vector<double> diff = util::VectorUtils<double>::ComputeDifference(q, current_point);
    double a = util::VectorUtils<double>::ComputeScalarProduct(derivative_dir1, derivative_dir1);
    double b = util::VectorUtils<double>::ComputeScalarProduct(derivative_dir2, derivative_dir1);
    double c = util::VectorUtils<double>::ComputeScalarProduct(diff, derivative_dir1);
    double d = util::VectorUtils<double>::ComputeScalarProduct(derivative_dir2, derivative_dir2);
    double e = util::VectorUtils<double>::ComputeScalarProduct(diff, derivative_dir2);
    double delta2 = (a * e - b * c) / (a * d - b * b);
    double delta1 = (c - b * delta2) / a;
    return std::array<double, 2>{delta1, delta2};
  }
};
}  // namespace spl

#endif  // SRC_SPL_PROJECTION_H_
