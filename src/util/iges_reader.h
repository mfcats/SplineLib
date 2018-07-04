/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_IGES_READER_H
#define SPLINELIB_IGES_READER_H

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "b_spline.h"
#include "spline.h"

namespace util {
class IGESReader {
 public:
  explicit IGESReader(std::string filename) : filename_(std::move(filename)) {}

  std::shared_ptr<spl::Spline<1>> ReadIGESFile() {
    std::ifstream newFile;
    newFile.open(filename_);
    std::string line;
    std::string globalSection;
    std::string parameterDataSection;
    while (getline(newFile, line)) {
      if (line[72] == 'G') {
        globalSection += line.substr(0, 72);
      }
      if (line[72] == 'P') {
        parameterDataSection += line.substr(0, 64);
      }
    }
    std::string delimiter = GetDelimiter(globalSection);
    std::string recordDelimiter = GetRecordDelimiter(globalSection);
    std::vector<std::string> data = SplitString(parameterDataSection, delimiter, recordDelimiter);
    if (data[0] == "406") {
      return Get1DSpline(data);
    }
    if (data[0] == "128") {
      //return Get2DSpline(data);
    }
    std::cout << std::endl << data[0] << std::endl;
    throw std::runtime_error( // NOLINT
        "The given IGES-file doesn't contain a readable B-Spline or NURBS.");
  }

 private:
  std::shared_ptr<spl::Spline<1>> Get1DSpline(const std::vector<std::string> &parameterData) {
    if (GetInteger(parameterData[5])) {
      return Get1DBSpline(parameterData);
    }
    //return Get1DNurbs(parameterData);
  }

/*  std::shared_ptr<spl::Spline> Get2DSpline(const std::vector<std::string> &parameterData) {
    if (GetInteger(parameterData[7])) {
      return Get2DBSpline(parameterData);
    }
    return Get2DNurbs(parameterData);
  }*/

  std::shared_ptr<spl::BSpline<1>> Get1DBSpline(const std::vector<std::string> &parameterData) {
    int upperSumIndex = GetInteger(parameterData[1]);
    std::array<int, 1> degree;
    degree[0] = GetInteger(parameterData[2]);
    std::vector<double> knots = GetData(7, upperSumIndex + degree[0] + 8, parameterData);
    std::vector<ParamCoord> knots_param;
    for(int i; i < knots.size(); ++i) {
      knots_param.push_back(ParamCoord{knots[i]});
    }
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector;
    knot_vector[0] = std::make_shared<baf::KnotVector>(knots_param);
    std::vector<baf::ControlPoint> control_points;
    for (int controlPoint = 0; controlPoint < upperSumIndex + 1; controlPoint++) {
      control_points.emplace_back(baf::ControlPoint(GetData(
          2 * upperSumIndex + degree[0] + 10 + 3 * controlPoint,
          2 * upperSumIndex + degree[0] + 10 + 3 * controlPoint + 2,
          parameterData)));
    }
    return std::make_shared<spl::BSpline<1>>(spl::BSplineGenerator<1>(knot_vector, degree, control_points));
  }

  /*std::shared_ptr<BSpline2D> Get2DBSpline(const std::vector<std::string> &parameterData) {
    int upperSumIndex1 = GetInteger(parameterData[1]);
    int upperSumIndex2 = GetInteger(parameterData[2]);
    int degree1 = GetInteger(parameterData[3]);
    int degree2 = GetInteger(parameterData[4]);
    KnotVector knotVector1 = KnotVector(GetData(10, upperSumIndex1 + degree1 + 11, parameterData));
    KnotVector knotVector2 = KnotVector(GetData(upperSumIndex1 + degree1 + 12,
                                                upperSumIndex1 + upperSumIndex2 + degree1 + degree2
                                                    + 13,
                                                parameterData));
    std::vector<ControlPoint> controlPoints;
    for (int controlPoint = 0; controlPoint < (upperSumIndex1 + 1) * (upperSumIndex2 + 1);
         controlPoint++) {
      controlPoints.emplace_back(ControlPoint(GetData(
          (upperSumIndex1 + 2) * (upperSumIndex2 + 2) + degree1 + degree2 + 11 + 3 * controlPoint,
          (upperSumIndex1 + 2) * (upperSumIndex2 + 2) + degree1 + degree2 + 11 + 3 * controlPoint + 2,
          parameterData)));
    }
    return std::make_shared<BSpline2D>(controlPoints,
                                       BasisFunctionDerivativeSet(knotVector1, degree1, degree1 + 1),
                                       BasisFunctionDerivativeSet(knotVector2, degree2, degree2 + 1));
  }

  std::shared_ptr<NURBS> Get1DNurbs(const std::vector<std::string> &parameterData) {
    int upperSumIndex = GetInteger(parameterData[1]);
    int degree = GetInteger(parameterData[2]);
    KnotVector knotVector = KnotVector(GetData(7, upperSumIndex + degree + 8, parameterData));
    std::vector<double>
        weights = GetData(upperSumIndex + degree + 9, 2 * upperSumIndex + degree + 9, parameterData);
    std::vector<WeightedControlPoint> controlPoints;
    for (int controlPoint = 0; controlPoint < upperSumIndex + 1; controlPoint++) {
      controlPoints.emplace_back(WeightedControlPoint(GetData(
          2 * upperSumIndex + degree + 10 + 3 * controlPoint,
          2 * upperSumIndex + degree + 10 + 3 * controlPoint + 2,
          parameterData), weights[controlPoint]));
    }
    return std::make_shared<NURBS>(controlPoints,
                                   BasisFunctionDerivativeSet(knotVector, degree, degree + 1));
  }

  std::shared_ptr<NURBS2D> Get2DNurbs(const std::vector<std::string> &parameterData) {
    int upperSumIndex1 = GetInteger(parameterData[1]);
    int upperSumIndex2 = GetInteger(parameterData[2]);
    int degree1 = GetInteger(parameterData[3]);
    int degree2 = GetInteger(parameterData[4]);
    KnotVector knotVector1 = KnotVector(GetData(10, upperSumIndex1 + degree1 + 11, parameterData));
    KnotVector knotVector2 = KnotVector(GetData(upperSumIndex1 + degree1 + 12,
                                                upperSumIndex1 + upperSumIndex2 + degree1 + degree2
                                                    + 13,
                                                parameterData));
    std::vector<double>
        weights = GetData(upperSumIndex1 + upperSumIndex2 + degree1 + degree2 + 14,
                          (upperSumIndex1 + 2) * (upperSumIndex2 + 2) + degree1 + degree2 + 10,
                          parameterData);
    std::vector<WeightedControlPoint> controlPoints;
    for (int controlPoint = 0; controlPoint < (upperSumIndex1 + 1) * (upperSumIndex2 + 1);
         controlPoint++) {
      controlPoints.emplace_back(WeightedControlPoint(GetData(
          (upperSumIndex1 + 2) * (upperSumIndex2 + 2) + degree1 + degree2 + 11 + 3 * controlPoint,
          (upperSumIndex1 + 2) * (upperSumIndex2 + 2) + degree1 + degree2 + 11 + 3 * controlPoint + 2,
          parameterData), weights[controlPoint]));
    }
    return std::make_shared<NURBS2D>(controlPoints,
                                     BasisFunctionDerivativeSet(knotVector1, degree1, degree1 + 1),
                                     BasisFunctionDerivativeSet(knotVector2, degree2, degree2 + 1));
  }*/

  std::vector<double> GetData(int start,
                                          int end,
                                          const std::vector<std::string> &parameterData) {
    std::vector<double> knots;
    for (int knot = start; knot <= end; knot++) {
      knots.emplace_back(GetDouble(parameterData[knot]));
    }
    return knots;
  }

  std::string GetDelimiter(const std::string &globalSection) {
    std::size_t found = globalSection.find_first_of('H');
    unsigned long number;
    std::istringstream(globalSection.substr(0, found)) >> number;
    if (number != 1) {
      throw std::runtime_error( // NOLINT
          "IGES-Files with a delimiter longer than one character can't be read.");
    }
    return globalSection.substr(found + 1, 1);
  }

  std::string GetRecordDelimiter(const std::string &globalSection) {
    std::size_t foundH = globalSection.find_first_of('H', globalSection.find_first_of('H') + 1);
    unsigned long number;
    std::istringstream(globalSection.substr(4, foundH - 4)) >> number;
    return globalSection.substr(foundH + 1, number);
  }

  int GetInteger(const std::string &string) {
    int number = 0;
    std::istringstream(string) >> number;
    return number;
  }

  double GetDouble(const std::string &string) {
    double number = 0;
    std::istringstream(string) >> number;
    return number;
  }

  std::vector<std::string> SplitString(std::string string,
                                                   const std::string &delimiter,
                                                   const std::string &recordDelimiter) {
    std::vector<std::string> splittedString;
    while (!string.empty()) {
      std::size_t foundDelimiter = string.find_first_of(delimiter);
      if (foundDelimiter == std::string::npos) {
        splittedString.push_back(string.substr(0, splittedString.back().find(recordDelimiter)));
        string.clear();
      } else {
        splittedString.push_back(string.substr(0, foundDelimiter));
        string.erase(0, foundDelimiter + 1);
      }
    }
    return splittedString;
  }

  std::string filename_;
};
}

#endif //SPLINELIB_IGES_READER_H
