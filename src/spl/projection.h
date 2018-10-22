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
#include "knot_vector.h"
#include "vector_utils.h"

namespace spl {
template<int DIM>
class Projection {
 public:
  static std::array<double, DIM> ProjectionOnSurface(const std::vector<double> &pointPhysicalCoords,
                                                     std::shared_ptr<spl::Spline<DIM>> spline,
                                                     std::array<int, 2> scattering = {10, 10}) {
    double tolerance = 0.00001;
    int iteration = 0;
    std::vector<int> dimensions;
    for (auto i = 0u; i < pointPhysicalCoords.size(); ++i) {
      dimensions.emplace_back(i);
    }
    std::array<ParamCoord, DIM>
        currentParamCoordGuess = FindInitialValue2D(scattering, pointPhysicalCoords, spline, dimensions);
    bool converged = false;

    while (!converged) {
      std::vector<double> derivative_dir1 = spline->EvaluateDerivative(currentParamCoordGuess, {0, 1, 2}, {1, 0});
      std::vector<double> derivative_dir2 = spline->EvaluateDerivative(currentParamCoordGuess, {0, 1, 2}, {0, 1});
      std::vector<double> currentPoint = spline->Evaluate(currentParamCoordGuess, {0, 1, 2});
      std::vector<double> q = GetQ(pointPhysicalCoords, currentPoint, derivative_dir1, derivative_dir2);
      std::array<double, 2> difference = GetDifferences(currentPoint, q, derivative_dir1, derivative_dir2);
      currentParamCoordGuess[0] = currentParamCoordGuess[0] + ParamCoord{difference[0]};
      currentParamCoordGuess[1] = currentParamCoordGuess[1] + ParamCoord{difference[1]};
      CheckParametricCoordinates(currentParamCoordGuess, spline);
      if (std::abs(difference[0]) < tolerance && std::abs(difference[1]) < tolerance) {
        converged = true;
      }
      ++iteration;
      if (iteration > 1000) {
        break;
      }
    }
    return {currentParamCoordGuess[0].get(), currentParamCoordGuess[1].get()};
  }

  static std::array<double, DIM> ProjectionOnSpline(const std::vector<double> &pointPhysicalCoords,
                                                    std::shared_ptr<spl::Spline<DIM>> spline) {
    double tolerance = 0.0001;
    int iteration = 0;
    std::vector<int> dimensions;
    for (auto i = 0u; i < pointPhysicalCoords.size(); ++i) {
      dimensions.emplace_back(i);
    }
    std::array<ParamCoord, DIM> currentParamCoordGuess = FindInitialValue(pointPhysicalCoords, spline, dimensions);
    bool converged = false;

    while (!converged) {
      if (DIM == 1) {
        std::vector<double> firstDer = spline->EvaluateDerivative(currentParamCoordGuess, dimensions, {1});
        std::vector<double> projectionVector =
            util::VectorUtils<double>::ComputeDifference(pointPhysicalCoords,
                                                         spline->Evaluate(currentParamCoordGuess, dimensions));

        // This is only the first order algorithm. An implemented but non-working version of the second order algorithm
        // can be found in commit 2ed993e6dcef3d184b70640f6b9498efae52747a.
        double delta = util::VectorUtils<double>::ComputeScalarProduct(firstDer, projectionVector)
            / util::VectorUtils<double>::ComputeScalarProduct(firstDer, firstDer);

        currentParamCoordGuess[0] = currentParamCoordGuess[0] + ParamCoord{delta};
        CheckParametricCoordinates(currentParamCoordGuess, spline);
        if (std::abs(delta) < tolerance) {
          converged = true;
        }
      }
      ++iteration;
      if (iteration > 1000) {
        break;
      }
    }
    return {currentParamCoordGuess[0].get()};
  }

 private:
  static std::array<ParamCoord, DIM> FindInitialValue(const std::vector<double> &pointPhysicalCoords,
                                                      const std::shared_ptr<spl::Spline<DIM>> &spline,
                                                      const std::vector<int> &dimensions) {
    std::vector<elm::Element> elements = spline->GetElementList();
    std::array<ParamCoord, DIM> paramCoords = {ParamCoord{0}};
    std::vector<double> splinePhysicalCoords =
        spline->Evaluate({ParamCoord{(0.5 * (elements[0].GetNode(1) - elements[0].GetNode(0)).get())}}, dimensions);
    double distance = util::VectorUtils<double>::ComputeTwoNorm(util::VectorUtils<double>::ComputeDifference(
        pointPhysicalCoords,
        splinePhysicalCoords));
    paramCoords[0] = ParamCoord{{0.5 * (elements[0].GetNode(1) - elements[0].GetNode(0)).get()}};
    for (auto i = 1u; i < elements.size(); ++i) {
      splinePhysicalCoords = spline->Evaluate({ParamCoord{
          0.5 * (elements[i].GetNode(1) - elements[i].GetNode(0)).get() + elements[i].GetNode(0).get()}}, dimensions);
      if (util::VectorUtils<double>::ComputeTwoNorm(util::VectorUtils<double>::ComputeDifference(pointPhysicalCoords,
                                                                                                 splinePhysicalCoords))
          < distance) {
        distance = util::VectorUtils<double>::ComputeTwoNorm(util::VectorUtils<double>::ComputeDifference(
            pointPhysicalCoords,
            splinePhysicalCoords));
        paramCoords[0] =
            ParamCoord{0.5 * (elements[i].GetNode(1) - elements[i].GetNode(0)).get() + elements[i].GetNode(0).get()};
      }
    }
    return paramCoords;
  }

  static void CheckParametricCoordinates(std::array<ParamCoord, DIM> &param_coords,
                                         const std::shared_ptr<spl::Spline<DIM>> &spline) {
    for (auto &coord : param_coords) {
      if (coord < spline->GetKnotVector(0)->GetKnot(0)) {
        coord = spline->GetKnotVector(0)->GetKnot(0);
      } else if (coord > spline->GetKnotVector(0)->GetLastKnot()) {
        coord = spline->GetKnotVector(0)->GetLastKnot();
      }
    }
  }

  static std::array<ParamCoord, 2> FindInitialValue2D(std::array<int, 2> scattering,
                                                      const std::vector<double> &pointPhysicalCoords,
                                                      const std::shared_ptr<spl::Spline<2>> &spline,
                                                      const std::vector<int> &dimensions) {
    double first_knot1 = spline->GetKnotVector(0)->GetKnot(0).get();
    double last_knot1 = spline->GetKnotVector(0)->GetLastKnot().get();
    double first_knot2 = spline->GetKnotVector(1)->GetKnot(0).get();
    double last_knot2 = spline->GetKnotVector(1)->GetLastKnot().get();
    std::array<ParamCoord, 2> param_coords = {ParamCoord(first_knot1), ParamCoord(first_knot2)};
    std::vector<double> physical_coords =
        spline->Evaluate({ParamCoord(first_knot1), ParamCoord(first_knot2)}, dimensions);
    double distance = util::VectorUtils<double>::ComputeTwoNorm(util::VectorUtils<double>::ComputeDifference(
        pointPhysicalCoords,
        physical_coords));
    for (int i = 1; i <= scattering[0]; ++i) {
      for (int j = 1; j <= scattering[1]; ++j) {
        ParamCoord coord1 = ParamCoord(i * (last_knot1 - first_knot1) / scattering[0]);
        ParamCoord coord2 = ParamCoord(j * (last_knot2 - first_knot2) / scattering[1]);
        physical_coords = spline->Evaluate({coord1, coord2}, dimensions);
        double current_dist = util::VectorUtils<double>::ComputeTwoNorm(util::VectorUtils<double>::ComputeDifference(
            pointPhysicalCoords,
            physical_coords));
        if (current_dist < distance) {
          param_coords = {coord1, coord2};
          distance = current_dist;
        }
      }
    }
    return param_coords;
  }

  static std::vector<double> GetQ(const std::vector<double> &projectionPoint, const std::vector<double> &currentPoint,
                                  const std::vector<double> &direction1, const std::vector<double> &direction2) {
    std::vector<double> normal_vector = util::VectorUtils<double>::CrossProduct(direction1, direction2);
    double length = util::VectorUtils<double>::ComputeTwoNorm(normal_vector);
    std::vector<double> normed_normal_vector = util::VectorUtils<double>::ScaleVector(normal_vector, 1 / length);
    double abstand = util::VectorUtils<double>::ComputeScalarProduct(currentPoint, normed_normal_vector)
        - util::VectorUtils<double>::ComputeScalarProduct(projectionPoint, normed_normal_vector);
    std::vector<double> scaled_normal_vector = util::VectorUtils<double>::ScaleVector(normed_normal_vector, -abstand);
    std::vector<double> q = util::VectorUtils<double>::ComputeDifference(projectionPoint, scaled_normal_vector);
    return q;
  }

  static std::array<double, 2> GetDifferences(const std::vector<double> &pointPhysicalCoords,
                                              const std::vector<double> &q,
                                              const std::vector<double> &derivativeInDirection1,
                                              const std::vector<double> &derivativeInDirection2) {
    std::vector<double> diff = util::VectorUtils<double>::ComputeDifference(q, pointPhysicalCoords);
    double a = util::VectorUtils<double>::ComputeScalarProduct(derivativeInDirection1, derivativeInDirection1);
    double b = util::VectorUtils<double>::ComputeScalarProduct(derivativeInDirection2, derivativeInDirection1);
    double c = util::VectorUtils<double>::ComputeScalarProduct(diff, derivativeInDirection1);
    double d = util::VectorUtils<double>::ComputeScalarProduct(derivativeInDirection2, derivativeInDirection2);
    double e = util::VectorUtils<double>::ComputeScalarProduct(diff, derivativeInDirection2);
    double delta2 = (a * e - b * c) / (a * d - b * b);
    double delta1 = (c - b * delta2) / a;
    return std::array<double, 2>{delta1, delta2};
  }
};
}  // namespace spl

#endif  // SRC_SPL_PROJECTION_H_
