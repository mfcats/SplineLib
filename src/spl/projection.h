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

namespace spl {
template<int DIM>
class Projection {
 public:
  Projection() {}

  std::vector<double> OrthogonalProjectionOntoCurve(std::vector<double> pointQ,
                                                    spl::Spline spline,
                                                    const std::vector<int> &dimensions) {
    std::vector<double> t = GetInitialValue(pointQ, spline);
    std::vector<double> dimensions;
    std::vector<double> derivative;
    std::vector<double> derivative2;
    for (int i = 0; i < DIM; ++i) {
      dimensions.emplace_back(i);
      derivative.emplace_back(1);
      derivative.emplace_back(2);
    }
    signDt = ComputeSignum(ComputeScalarProduct(spline.EvaluateDerivative(t, dimensions, derivative, pointQ),
                                                ComputePiecewiseVectorDifference(pointQ,
                                                                                 spline.Evaluate(t, dimensions))));
    while (dt >= dtEnd) {
      double dt = sqrt(2 * ComputeArea(spline.EvaluateDerivative(t, dimensions, derivative),
                                       ComputePiecewiseVectorDifference(pointQ, spline.Evaluate(t, dimensions)))
                           / ComputeArea(spline.EvaluateDerivative(t, dimensions, derivative),
                                         spline.EvaluateDerivative(t, dimensions, derivative2)));
      t += signDt * dt;
    }
    std::vector<double> pointP = spline.Evaluate(t, dimensions);
    std::transform(pointQ.begin(), pointQ.end(), pointP.begin(), pointP.begin(), std::minus<double>());
    std::transform(pointP.begin(), pointP.end(), pointP.begin(), pointP.begin(), std::multiplies<double>());
    return sqrt(std::accumulate(pointP.begin(), pointP.end(), 0));
  }

  int ComputeSignum(double x) {
    return ((x > 0) ? 1 : ((x < 0) ? -1 : 0));
  }

  std::vector<double> ComputePiecewiseVectorDifference(std::vector<double> vectorA, std::vector<double> vectorB) {
    std::transform(vectorA.begin(), vectorA.end(), vectorB.begin(), vectorB.begin(), std::minus<double>());
    return vectorB;
  }

  std::vector<double> ComputeArea(std::vector<double> vectorA, std::vector<double> vectorB) {
    std::vector<double> area;
    if (vectorA.Size() == 2) {
      area.emplace_back(vectorA[0] * vectorB[1] - vectorA[1] * vectorB[0]);
    } else if (vectorA.Size() == 3) {
      area.emplace_back(vectorA[1] * vectorB[2] - vectorA[2] * vectorB[1]);
      area.emplace_back(vectorA[2] * vectorB[0] - vectorA[0] * vectorB[2]);
      area.emplace_back(vectorA[0] * vectorB[1] - vectorA[1] * vectorB[0]);
    }
    return area;
  }

  double ComputeScalarProduct(std::vector<double> vectorA, std::vector<double> vectorB) {
    double sum = 0;
    for (int i = 0; i < vectorA.Size(); ++i) {
      double product = 0;
      product = vectorA[i] * vectorB[i];
      sum += product;
    }
    return sum;
  }
};
}

#endif //SPLINELIB_PROJECTION_H
