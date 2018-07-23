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
            /util::VectorUtils<double>::ComputeScalarProduct(firstDer, firstDer);

        projectionPointParamCoords[0] = projectionPointParamCoords[0] + ParamCoord{delta};
        if (projectionPointParamCoords[0] < spline->GetKnotVector(0).GetKnot(0)) {
          projectionPointParamCoords[0] = spline->GetKnotVector(0).GetKnot(0);
        } else if (projectionPointParamCoords[0] > spline->GetKnotVector(0).GetLastKnot()) {
          projectionPointParamCoords[0] = spline->GetKnotVector(0).GetLastKnot();
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

  static std::array<ParamCoord, DIM> FindInitialValue(const std::vector<double> &pointPhysicalCoords,
                                                      const std::unique_ptr<spl::Spline<DIM>> &spline,
                                                      const std::vector<int> &dimensions) {
    std::vector<elm::Element> elements = spline->GetElementList();
    std::array<ParamCoord, DIM> paramCoords = {ParamCoord{0}};
    std::vector<double> splinePhysicalCoords =
        spline->Evaluate({ParamCoord{(0.5*(elements[0].GetNode(1) - elements[0].GetNode(0)))}}, dimensions);
    double distance = util::VectorUtils<double>::ComputeTwoNorm(util::VectorUtils<double>::ComputeDifference(
        pointPhysicalCoords,
        splinePhysicalCoords));
    paramCoords[0] = ParamCoord{{0.5*(elements[0].GetNode(1) - elements[0].GetNode(0))}};
    for (auto i = 1u; i < elements.size(); ++i) {
      splinePhysicalCoords = spline->Evaluate({ParamCoord{
          0.5*(elements[i].GetNode(1) - elements[i].GetNode(0)) + elements[i].GetNode(0)}}, dimensions);
      if (util::VectorUtils<double>::ComputeTwoNorm(util::VectorUtils<double>::ComputeDifference(pointPhysicalCoords,
                                                                                                 splinePhysicalCoords))
          < distance) {
        distance = util::VectorUtils<double>::ComputeTwoNorm(util::VectorUtils<double>::ComputeDifference(
            pointPhysicalCoords,
            splinePhysicalCoords));
        paramCoords[0] = ParamCoord{0.5*(elements[i].GetNode(1) - elements[i].GetNode(0)) + elements[i].GetNode(0)};
      }
    }
    return paramCoords;
  }
};
}  // namespace spl

#endif  // SRC_SPL_PROJECTION_H_
