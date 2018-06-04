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
#include <iostream>

#include "b_spline.h"
#include "numeric_settings.h"

namespace spl {
class Projection {
 public:
  static std::array<double, 1> ProjectionOnSpline(std::vector<double> pointPhysicalCoords,
                                              spl::BSpline<1> *spline) {
    double kappa;
    double tolerance = 0.0001;
    int iteration = 0;
    double distance;
    double delta;
    int signum;
    std::array<double, 1> projectionPointParamCoords = {0.5};  //u  //array length DIM
    std::vector<int> dimensions;
    for (int i = 0; i < pointPhysicalCoords.size(); ++i) {
      dimensions.emplace_back(i);
    }
    bool converged = false;

    while (not converged) {
      if (true /*spline->parameter_space_.size() == 1*/) {   //spline dimension is one
        std::vector<double> firstDer = spline->EvaluateDerivative(projectionPointParamCoords, dimensions, {1});
        std::vector<double> secondDer = spline->EvaluateDerivative(projectionPointParamCoords, dimensions, {2});
        kappa = ComputeArea(firstDer, secondDer) / pow(ComputeTwoNorm(firstDer), 3);
        std::vector<double> projectionVector = ComputePiecewiseVectorDifference(pointPhysicalCoords, spline->Evaluate(projectionPointParamCoords, dimensions));

        if (/*kappa >= 10e-8*/ false) {   //second order algorithm
          delta = sqrt(2 * ComputeArea(firstDer, projectionVector) / (kappa * pow(ComputeTwoNorm(firstDer), 3)));
          signum = ComputeScalarProduct(firstDer, projectionVector) / abs(ComputeScalarProduct(firstDer, projectionVector));
        } else {
          delta = ComputeScalarProduct(firstDer, projectionVector) / ComputeScalarProduct(firstDer, firstDer);
          signum = 1;
        }
        projectionPointParamCoords[0] += signum * delta;
        if (projectionPointParamCoords[0] < spline->GetKnotVector(0).knot(0)) {
          projectionPointParamCoords[0] = spline->GetKnotVector(0).knot(0);
        } else if (projectionPointParamCoords[0] > spline->GetKnotVector(0).GetLastKnot()) {
          projectionPointParamCoords[0] = spline->GetKnotVector(0).GetLastKnot();
        }
        if (abs(delta) < tolerance) {
          converged = true;
        }
      }
      std::cout << signum * delta << "   " << spline->Evaluate(projectionPointParamCoords, dimensions)[0] << "   " << spline->Evaluate(projectionPointParamCoords, dimensions)[1] << std::endl;
      ++iteration;
      if (iteration > 1000) {
        break;
      }
    }
    return projectionPointParamCoords;
  }

  static int ComputeSignum(double x) {
    return ((x > 0) ? 1 : ((x < 0) ? -1 : 0));
  }

  static double ComputeTwoNorm(std::vector<double> vectorA) {
    std::transform(vectorA.begin(), vectorA.end(), vectorA.begin(), vectorA.begin(), std::multiplies<double>());
    return sqrt(std::accumulate(vectorA.begin(), vectorA.end(), 0));
  }

  static std::vector<double> ComputePiecewiseVectorDifference(std::vector<double> vectorA, std::vector<double> vectorB) {
    std::transform(vectorA.begin(), vectorA.end(), vectorB.begin(), vectorB.begin(), std::minus<double>());
    return vectorB;
  }

  static double ComputeArea(std::vector<double> vectorA, std::vector<double> vectorB) {
    return sqrt(abs(ComputeScalarProduct(vectorA, vectorA) *
        ComputeScalarProduct(vectorB, vectorB) - pow(ComputeScalarProduct(vectorA, vectorB), 2)));
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
