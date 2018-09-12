/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_IGES_2D_BSPLINE_GENERATOR_H_
#define SRC_SPL_IGES_2D_BSPLINE_GENERATOR_H_

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
class IGES2DBSplineGenerator : public BSplineGenerator<2> {
 public:
  IGES2DBSplineGenerator() {}

  void ReadIGESFile(const char* filename, int entityToBeRead) {
    std::ifstream newFile;
    newFile.open(filename);
    if (!newFile.good()) {
      throw std::runtime_error("IGES file could not be opened.");
    }
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
    if ((parameterData[0] == 128) && (parameterData[7] == 1)) {
      std::array<int, 2> upperSumIndex;
      upperSumIndex[0] = static_cast<int>(parameterData[1]);
      upperSumIndex[1] = static_cast<int>(parameterData[2]);
      std::array<Degree, 2> degree;
      degree[0] = Degree{static_cast<int>(parameterData[3])};
      degree[1] = Degree{static_cast<int>(parameterData[4])};
      std::array<std::shared_ptr<baf::KnotVector>, 2> knot_vector;
      std::vector<double> weights;
      std::vector<baf::ControlPoint> control_points;
      std::array<std::array<int, 2>, 2> knotsStartEnd;
      std::array<int, 2> weightsStartEnd;
      std::array<int, 2> controlPointsStartEnd;
      knotsStartEnd[0][0] = 10;
      knotsStartEnd[0][1] = knotsStartEnd[0][0] + upperSumIndex[0] + degree[0].get() + 1;
      knotsStartEnd[1][0] = knotsStartEnd[0][1] + 1;
      knotsStartEnd[1][1] = knotsStartEnd[1][0] + upperSumIndex[1] + degree[1].get() + 1;
      weightsStartEnd[0] = knotsStartEnd[1][1] + 1;
      weightsStartEnd[1] = weightsStartEnd[0] - 1 + ((1 + upperSumIndex[0]) * (1 + upperSumIndex[1]));
      controlPointsStartEnd[0] = weightsStartEnd[1] + 1;
      controlPointsStartEnd[1] = controlPointsStartEnd[0] - 3 + (3 * (1 + upperSumIndex[0]) * (1 + upperSumIndex[1]));
      std::array<std::vector<ParamCoord>, 2> knots;
      for (int i = knotsStartEnd[0][0]; i <= knotsStartEnd[0][1]; ++i) {
        knots[0].push_back(ParamCoord{parameterData[i]});
      }
      for (int i = knotsStartEnd[1][0]; i <= knotsStartEnd[1][1]; ++i) {
        knots[1].push_back(ParamCoord{parameterData[i]});
      }
      knot_vector[0] = std::make_shared<baf::KnotVector>(knots[0]);
      knot_vector[1] = std::make_shared<baf::KnotVector>(knots[1]);
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
      std::array<int, 2> number_of_points;
      for (int i = 0; i < 2; ++i) {
        number_of_points[i] = knot_vector[i]->GetNumberOfKnots() - degree[i].get() - 1;
      }
      this->physical_space_ = std::make_shared<PhysicalSpace<2>>(control_points, number_of_points);
      this->parameter_space_ = std::make_shared<ParameterSpace<2>>(knot_vector, degree);
    } else {
      throw std::runtime_error("You are trying to read an entity of the wrong type.");
    }
  }

  std::array<int, 2> GetParameterSectionStartEndPointers(std::vector<std::string> directoryEntrySection,
                                                         int entityToBeRead) {
    std::string parameterDataStartPointer = trim(directoryEntrySection[entityToBeRead * 2].substr(8, 8));
    std::string parameterDataLineCount = trim(directoryEntrySection[entityToBeRead * 2 + 1].substr(24, 8));
    std::array<int, 2> ParameterSectionStartEndPointers;
    ParameterSectionStartEndPointers[0] = GetInteger(trim(parameterDataStartPointer));
    ParameterSectionStartEndPointers[1] =
        GetInteger(trim(parameterDataStartPointer)) + GetInteger(trim(parameterDataLineCount)) - 1;
    return ParameterSectionStartEndPointers;
  }

  std::vector<double> ParameterSectionToVector(std::vector<std::string> parameterSection,
                                               std::array<int, 2> ParameterSectionStartEndPointers) {
    std::vector<double> parameterSectionVector;
    int first = ParameterSectionStartEndPointers[0] - 1;
    int last = ParameterSectionStartEndPointers[1] - 1;
    for (int i = first; i <= last; ++i) {
      auto temp = DelimitedStringToVector(parameterSection[i]);
      for (int j = 0; j < temp.size(); ++j) {
        parameterSectionVector.push_back(temp[j]);
      }
    }
    return parameterSectionVector;
  }

  std::vector<double> DelimitedStringToVector(std::string str) {
    std::vector<double> vector;
    std::size_t found1;
    std::size_t found2;
    while (!str.empty()) {
      found1 = str.find_first_of(',');
      found2 = str.find_first_of(';');
      if ((found1 < found2) && (found1 != 0)) {
        vector.push_back(GetDouble(str.substr(0, found1)));
        str.erase(0, found1 + 1);
      } else if ((found2 < found1) && (found2 != 0)) {
        vector.push_back(GetDouble(str.substr(0, found2)));
        str.erase(0, found2 + 1);
      } else {
        str.erase(0, 1);
      }
    }
    return vector;
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

  double GetDouble(const std::string &string) {
    double number = 0;
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
};
}  // namespace spl

#endif  // SRC_SPL_IGES_2D_BSPLINE_GENERATOR_H_
