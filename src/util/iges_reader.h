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

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "b_spline.h"
#include "nurbs.h"
#include "spline.h"

namespace util {
class IGESReader {
 public:
  explicit IGESReader(std::string filename) : filename_(std::move(filename)) {}

  std::shared_ptr<spl::Spline<1>> ReadIGESFile(int entityToBeRead) {
    std::ifstream newFile;
    newFile.open(filename_);
    std::string line;
    std::string globalSection;
    std::vector<std::string> directoryEntrySection;
    std::vector<std::string> parameterDataSection;
    while (getline(newFile, line)) {
      if (line[72] == 'D') {
        directoryEntrySection.push_back(line.substr(0, 72));
      }
      if (line[72] == 'P') {
        parameterDataSection.push_back(line.substr(0, 64));
      }
    }
    return GetSpline(ParameterSectionToVector(parameterDataSection, GetParameterSectionStartEndPointers(directoryEntrySection, entityToBeRead)));
  }

 private:
  std::shared_ptr<spl::Spline<1>> GetSpline(const std::vector<double> &parameterData) {
    if (parameterData[0] == 126) {
      int upperSumIndex = int(parameterData[1]);
      std::array<int, 1> degree;
      degree[0] = int(parameterData[2]);
      std::vector<ParamCoord> knots;
      for (int i = 7; i <= (upperSumIndex + degree[0] + 7); ++i) {
        knots.push_back(ParamCoord{parameterData[i]});
      }
      std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector;
      knot_vector[0] = std::make_shared<baf::KnotVector>(knots);
      std::vector<baf::ControlPoint> control_points;
      for (int i = 0; i <= upperSumIndex; ++i) {
        std::vector<double> controlPoint;
        for (int j = (2*upperSumIndex + degree[0] + 9 + 3*i); j <= (2*upperSumIndex + degree[0] + 9 + 3*i + 2); ++j) {
          controlPoint.push_back(parameterData[i]);
        }
        control_points.push_back(baf::ControlPoint(controlPoint));
      }
      if (parameterData[4]) {
        std::vector<double> weights;
        for (int i = (upperSumIndex + degree[0] + 8); i <= (2*upperSumIndex + degree[0] + 8); ++i) {
          weights.push_back(parameterData[i]);
          return std::make_shared<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
        }
      }
    } else {
      throw std::runtime_error("The given IGES-file doesn't contain a readable B-Spline or NURBS.");
    }
  }

  std::array<int, 2> GetParameterSectionStartEndPointers(std::vector<std::string> directoryEntrySection, int entityToBeRead) {
    std::string parameterDataStartPointer = directoryEntrySection[entityToBeRead*2];
    std::string parameterDataLineCount = directoryEntrySection[entityToBeRead*2 + 1];
    parameterDataStartPointer = trim(parameterDataStartPointer);
    parameterDataLineCount = trim(parameterDataLineCount);
    parameterDataStartPointer = parameterDataStartPointer.substr(8,8);
    parameterDataLineCount = parameterDataLineCount.substr(24,8);
    trim(parameterDataStartPointer);
    trim(parameterDataLineCount);
    std::array<int, 2> ParameterSectionStartEndPointers;
    ParameterSectionStartEndPointers[0] = GetInteger(parameterDataStartPointer);
    ParameterSectionStartEndPointers[1] = GetInteger(parameterDataStartPointer) + GetInteger(parameterDataLineCount) - 1;
    return ParameterSectionStartEndPointers;
  };

  std::vector<double> ParameterSectionToVector(std::vector<std::string> parameterSection, std::array<int, 2> ParameterSectionStartEndPointers) {
    std::vector<double> parameterSectionVector;
    int first = ParameterSectionStartEndPointers[0] - 1;
    int last = ParameterSectionStartEndPointers[1] - 1;
    for (int i = first; i <= last; ++i) {
      std::stringstream ss(parameterSection[i]);
      double d;
      while (ss >> d) {
        parameterSectionVector.push_back(d);
        if (ss.peek() == ',') {
          ss.ignore();
        }
      }
    }
    return parameterSectionVector;
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

  inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
  {
    s.erase(0, s.find_first_not_of(t));
    return s;
  }

  inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
  {
    s.erase(s.find_last_not_of(t) + 1);
    return s;
  }

  inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v")
  {
    return ltrim(rtrim(s, t), t);
  }

  std::string filename_;
};
}

#endif //SPLINELIB_IGES_READER_H
