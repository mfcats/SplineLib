/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_IGES_1D_BSPLINE_GENERATOR_H
#define SPLINELIB_IGES_1D_BSPLINE_GENERATOR_H

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
#include "spline.h"
#include "spline_generator.h"
#include "b_spline_generator.h"

namespace spl {
class IGES1DBSplineGenerator : public BSplineGenerator<1> {
 public:
  explicit IGES1DBSplineGenerator(std::string filename) : filename_(std::move(filename)) {}

  void ReadIGESFile(int entityToBeRead) {
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
    ReadParameterData(ParameterSectionToVector(parameterDataSection,
                                               GetParameterSectionStartEndPointers(directoryEntrySection,
                                                                                   entityToBeRead)));
  }

 private:
  void ReadParameterData(const std::vector<double> &parameterData) {
    if ((parameterData[0] == 126) && (parameterData[5] == 1)) {
      int upperSumIndex = int(parameterData[1]);
      std::array<int, 1> degree;
      degree[0] = int(parameterData[2]);
      std::array<baf::KnotVector, 1> knot_vector;
      std::vector<double> weights;
      std::vector<baf::ControlPoint> control_points;
      std::array<int, 2> knotsStartEnd;
      std::array<int, 2> weightsStartEnd;
      std::array<int, 2> controlPointsStartEnd;
      knotsStartEnd[0] = 7;
      knotsStartEnd[1] = knotsStartEnd[0] + upperSumIndex + degree[0] + 1;
      weightsStartEnd[0] = knotsStartEnd[1] + 1;
      weightsStartEnd[1] = weightsStartEnd[0] + upperSumIndex;
      controlPointsStartEnd[0] = weightsStartEnd[1] + 1;
      controlPointsStartEnd[1] = controlPointsStartEnd[0] + (3 * upperSumIndex);
      std::vector<ParamCoord> knots;
      for (int i = knotsStartEnd[0]; i <= knotsStartEnd[1]; ++i) {
        knots.push_back(ParamCoord{parameterData[i]});
      }
      knot_vector[0] = baf::KnotVector(knots);
      for (int i = weightsStartEnd[0]; i <= weightsStartEnd[1]; ++i) {
        weights.push_back(parameterData[i]);
      }
      std::vector<double> controlPointCoordinates;
      for (int i = controlPointsStartEnd[0]; i <= controlPointsStartEnd[1]; ++i) {
        controlPointCoordinates.push_back(parameterData[i]);
      }
      for (int i = 0; i < controlPointCoordinates.size(); i += 3) {
        control_points.push_back(baf::ControlPoint({controlPointCoordinates[i],
                                                    controlPointCoordinates[i + 1],
                                                    controlPointCoordinates[i + 2]}));
      }
      std::array<int, 1> number_of_points;
      for (int i = 0; i < 1; ++i) {
        number_of_points[i] = knot_vector[i].GetNumberOfKnots() - degree[i] - 1;
      }
      this->physical_space_ptr = std::make_shared<PhysicalSpace<1>>(control_points, number_of_points);
      this->parameter_space_ptr = std::make_shared<ParameterSpace<1>>(knot_vector, degree);
  } else {
      throw std::runtime_error("You are trying to read an entity of the wrong type.");
    }
  }

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

  static inline std::string trim(std::string s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
      return !std::isspace(ch);
    }));
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
      return !std::isspace(ch);
    }).base(), s.end());
    return s;
  }

  std::string filename_;
};
}

#endif //SPLINELIB_IGES_1D_BSPLINE_GENERATOR_H