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
    return Get1DSpline(ParameterSectionToVector(parameterDataSection, GetParameterSectionStartEndPointers(directoryEntrySection, entityToBeRead)));
  }

 private:
  std::shared_ptr<spl::Spline<1>> Get1DSpline(const std::vector<double> &parameterData) {
    if (parameterData[0] == 126) {
      if (parameterData[5] == 1) {
        return Get1DBSpline(parameterData);
      } else if (parameterData[5] == 0) {
        return Get1DNURBS(parameterData);
      }
    } else if (parameterData[0] == 128) {
      if (parameterData[7] == 1) {
        //return Get2DBSpline(parameterData);
      } else if (parameterData[7] == 0) {
        //return Get2DNURBS(parameterData);
      }
    } else {
      throw std::runtime_error("The given IGES-file doesn't contain a readable B-Spline or NURBS.");
    }
  }

  std::shared_ptr<spl::BSpline<1>> Get1DBSpline(const std::vector<double> &parameterData) {
    int upperSumIndex = int(parameterData[1]);
    std::array<int, 1> degree;
    degree[0] = int(parameterData[2]);
    std::vector<ParamCoord> knots;
    for (int i = 7; i <= (upperSumIndex + degree[0] + 8); ++i) {
      knots.push_back(ParamCoord{parameterData[i]});
    }
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector;
    knot_vector[0] = std::make_shared<baf::KnotVector>(knots);
    std::vector<baf::ControlPoint> control_points;
    for (int i = 0; i <= upperSumIndex; ++i) {
      control_points.emplace_back(baf::ControlPoint(ExtractPartOfVector(
          2 * upperSumIndex + degree[0] + 10 + 3 * i,
          2 * upperSumIndex + degree[0] + 10 + 3 * i + 2,
          parameterData)));
    }
      return std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
    }

  std::shared_ptr<spl::NURBS<1>> Get1DNURBS(const std::vector<double> &parameterData) {
    int upperSumIndex = int(parameterData[1]);
    std::array<int, 1> degree;
    degree[0] = int(parameterData[2]);
    std::vector<ParamCoord> knots;
    for (int i = 7; i <= (upperSumIndex + degree[0] + 8); ++i) {
      knots.push_back(ParamCoord{parameterData[i]});
    }
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector;
    knot_vector[0] = std::make_shared<baf::KnotVector>(knots);
    std::vector<baf::ControlPoint> control_points;
    for (int i = 0; i <= upperSumIndex; ++i) {
      control_points.emplace_back(baf::ControlPoint(ExtractPartOfVector(
          2 * upperSumIndex + degree[0] + 10 + 3 * i,
          2 * upperSumIndex + degree[0] + 10 + 3 * i + 2,
          parameterData)));
    std::vector<double> weights;
    for (int i = (upperSumIndex + degree[0] + 8); i <= (2*upperSumIndex + degree[0] + 8); ++i) {
      weights.push_back(parameterData[i]);
    }
    return std::make_shared<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
    }
  }

  //std::shared_ptr<spl::BSpline<2>> Get2DBSpline(const std::vector<double> &parameterData) {}

  //std::shared_ptr<spl::NURBS<2>> Get2DNURBS(const std::vector<double> &parameterData) {}

  std::array<int, 2> GetParameterSectionStartEndPointers(std::vector<std::string> directoryEntrySection, int entityToBeRead) {
    std::string parameterDataStartPointer = trim(directoryEntrySection[entityToBeRead*2]).substr(8,8);
    std::string parameterDataLineCount = trim(directoryEntrySection[entityToBeRead*2 + 1]).substr(24,8);
    std::array<int, 2> ParameterSectionStartEndPointers;
    ParameterSectionStartEndPointers[0] = GetInteger(trim(parameterDataStartPointer));
    ParameterSectionStartEndPointers[1] = GetInteger(trim(parameterDataStartPointer)) + GetInteger(trim(parameterDataLineCount)) - 1;
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

  std::vector<double> ExtractPartOfVector(int start,
                              int end,
                              const std::vector<double> &input) {
    std::vector<double> output;
    for (int i = start; i <= end; ++i) {
      output.emplace_back(input[i]);
    }
    return output;
  }

  int GetInteger(const std::string &string) {
    int number = 0;
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
