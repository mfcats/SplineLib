/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_PROJECTION_H
#define SPLINELIB_PROJECTION_H

#include <vector>

#include "b_spline.h"
#include "numeric_settings.h"

namespace spl {
class Projection {
 public:
  static double OrthogonalProjectionOntoCurve(std::vector<double> pointQ,
                                              std::unique_ptr<spl::BSpline<1>> spline) {
    double t = 0.5;//GetInitialValue(pointQ, spline);
    std::vector<int> dimensions;
    for (int i = 0; i < pointQ.size(); ++i) {
      dimensions.emplace_back(i);
    }
    int signDt = Projection::ComputeSignum(Projection::ComputeScalarProduct(spline->EvaluateDerivative({t}, dimensions, {1}),
                                                Projection::ComputePiecewiseVectorDifference(pointQ,
                                                                                 spline->Evaluate({t}, dimensions))));
    double dtEnd = util::NumericSettings<double>::kEpsilon();
    double dt = dtEnd;
    while (dt >= dtEnd) {
      dt = sqrt(2 * Projection::ComputeArea(spline->EvaluateDerivative({t}, dimensions, {1}),
                                       Projection::ComputePiecewiseVectorDifference(pointQ, spline->Evaluate({t}, dimensions)))
                           / Projection::ComputeArea(spline->EvaluateDerivative({t}, dimensions, {1}),
                                         spline->EvaluateDerivative({t}, dimensions, {2})));
      t += signDt * dt;
    }
    //std::vector<double> pointP = spline.Evaluate(t, dimensions);
    //std::transform(pointQ.begin(), pointQ.end(), pointP.begin(), pointP.begin(), std::minus<double>());
    //std::transform(pointP.begin(), pointP.end(), pointP.begin(), pointP.begin(), std::multiplies<double>());
    //return sqrt(std::accumulate(pointP.begin(), pointP.end(), 0));
    return t;
  }

  static int ComputeSignum(double x) {
    return ((x > 0) ? 1 : ((x < 0) ? -1 : 0));
  }

  static std::vector<double> ComputePiecewiseVectorDifference(std::vector<double> vectorA, std::vector<double> vectorB) {
    std::transform(vectorA.begin(), vectorA.end(), vectorB.begin(), vectorB.begin(), std::minus<double>());
    return vectorB;
  }

  static double ComputeArea(std::vector<double> vectorA, std::vector<double> vectorB) {
    double area;
    if (vectorA.size() == 2) {
      area = vectorA[0] * vectorB[1] - vectorA[1] * vectorB[0];
    } else if (vectorA.size() == 3) {
      area = (vectorA[1] * vectorB[2] - vectorA[2] * vectorB[1]) + (vectorA[2] * vectorB[0] - vectorA[0] * vectorB[2]) + (vectorA[0] * vectorB[1] - vectorA[1] * vectorB[0]);
    }
    return area;
  }

  static double ComputeScalarProduct(std::vector<double> vectorA, std::vector<double> vectorB) {
    double sum = 0;
    for (int i = 0; i < vectorA.size(); ++i) {
      double product = 0;
      product = vectorA[i] * vectorB[i];
      sum += product;
    }
    return sum;
  }
};
}

#endif //SPLINELIB_PROJECTION_H
