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
                                                     std::unique_ptr<spl::Spline<DIM>> spline) {
    std::array<ParamCoord, DIM> currentParamCoordGuess = {ParamCoord(0.6), ParamCoord(0.9)};
    double tolerance = 0.00001;
    int iteration = 0;
    bool converged = false;

    while (!converged) {
      std::vector<double>
          derivativeInDirection1 = spline->EvaluateDerivative(currentParamCoordGuess, {0, 1, 2}, {1, 0});
      std::vector<double>
          derivativeInDirection2 = spline->EvaluateDerivative(currentParamCoordGuess, {0, 1, 2}, {0, 1});
      std::vector<double>
          normal_vector = util::VectorUtils<double>::CrossProduct(derivativeInDirection1, derivativeInDirection2);
      std::vector<double> currentPoint = spline->Evaluate(currentParamCoordGuess, {0, 1, 2});
      std::vector<double> q = GetQ(pointPhysicalCoords, currentPoint, derivativeInDirection1, derivativeInDirection2);
      auto a = q[0];
      auto b = q[1];
      auto c = q[2];
      std::array<double, 2>
          difference = GetDifferences(currentPoint, q, derivativeInDirection1, derivativeInDirection2);
      currentParamCoordGuess[0] = currentParamCoordGuess[0] + ParamCoord{difference[0]};
      currentParamCoordGuess[1] = currentParamCoordGuess[1] + ParamCoord{difference[1]};
      if (std::abs(difference[0]) < tolerance && std::abs(difference[1]) < tolerance) {
        converged = true;
      }
      ++iteration;
      if (iteration > 1000) {
        break;
      }
    }
    auto r = spline->Evaluate(currentParamCoordGuess, {0})[0];
    auto s = spline->Evaluate(currentParamCoordGuess, {1})[0];
    auto t = spline->Evaluate(currentParamCoordGuess, {2})[0];

    auto u = spline->Evaluate({ParamCoord(0.8614466), ParamCoord(0.556805)}, {0})[0];
    auto v = spline->Evaluate({ParamCoord(0.8614466), ParamCoord(0.556805)}, {1})[0];
    auto w = spline->Evaluate({ParamCoord(0.8614466), ParamCoord(0.556805)}, {2})[0];
    return {currentParamCoordGuess[0].get(), currentParamCoordGuess[1].get()};
  }

  static std::array<double, DIM> ProjectionOnSpline(const std::vector<double> &pointPhysicalCoords,
                                                    std::unique_ptr<spl::Spline<DIM>> spline) {
    double tolerance = 0.0001;
    int iteration = 0;
    std::vector<int> dimensions;
    for (auto i = 0u; i < pointPhysicalCoords.size(); ++i) {
      dimensions.emplace_back(i);
    }
    std::array<ParamCoord, DIM> projectionPointParamCoords = FindInitialValue(pointPhysicalCoords, spline, dimensions);
    bool converged = false;

    while (!converged) {
      if (DIM == 1) {
        std::vector<double> firstDer = spline->EvaluateDerivative(projectionPointParamCoords, dimensions, {1});
        std::vector<double> projectionVector =
            util::VectorUtils<double>::ComputeDifference(pointPhysicalCoords,
                                                         spline->Evaluate(projectionPointParamCoords, dimensions));

        // This is only the first order algorithm. An implemented but non-working version of the second order algorithm
        // can be found in commit 2ed993e6dcef3d184b70640f6b9498efae52747a.
        double delta = util::VectorUtils<double>::ComputeScalarProduct(firstDer, projectionVector)
            / util::VectorUtils<double>::ComputeScalarProduct(firstDer, firstDer);

        projectionPointParamCoords[0] = projectionPointParamCoords[0] + ParamCoord{delta};
        if (projectionPointParamCoords[0].get() < spline->GetKnotVector(0)->GetKnot(0).get()) {
          projectionPointParamCoords[0] = spline->GetKnotVector(0)->GetKnot(0);
        } else if (projectionPointParamCoords[0] > spline->GetKnotVector(0)->GetLastKnot()) {
          projectionPointParamCoords[0] = spline->GetKnotVector(0)->GetLastKnot();
        }
        if (std::abs(delta) < tolerance) {
          converged = true;
        }
      }
      ++iteration;
      if (iteration > 1000) {
        break;
      }
    }
    return {projectionPointParamCoords[0].get()};
  }

 private:
  static std::array<ParamCoord, DIM> FindInitialValue(const std::vector<double> &pointPhysicalCoords,
                                                      const std::unique_ptr<spl::Spline<DIM>> &spline,
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
    double x = diff[0];
    double y = diff[1];
    double z = diff[2];
    double a = util::VectorUtils<double>::ComputeScalarProduct(derivativeInDirection1, derivativeInDirection1);
    double b = util::VectorUtils<double>::ComputeScalarProduct(derivativeInDirection2, derivativeInDirection1);
    double c = util::VectorUtils<double>::ComputeScalarProduct(diff, derivativeInDirection1);
    double d = util::VectorUtils<double>::ComputeScalarProduct(derivativeInDirection1, derivativeInDirection2);
    double e = util::VectorUtils<double>::ComputeScalarProduct(derivativeInDirection2, derivativeInDirection2);
    double f = util::VectorUtils<double>::ComputeScalarProduct(diff, derivativeInDirection2);
    double delta2 = (f - d * c / a) / (e - d * b / a);
    double delta1 = (c - b * f / e) / (a - b * d / e);
    return std::array<double, 2>{delta1, delta2};
  }
};
}  // namespace spl

#endif  // SRC_SPL_PROJECTION_H_
