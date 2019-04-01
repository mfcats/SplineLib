/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "iges_reader.h"

#include <fstream>
#include <sstream>
#include <utility>

#include "b_spline.h"
#include "nurbs.h"
#include "string_operations.h"

std::vector<std::any> io::IGESReader::ReadFile(const char *filename) {
  std::ifstream newFile;
  newFile.open(filename);
  if (!newFile.good()) {
    throw std::runtime_error("IGES file could not be opened.");
  }
  std::string line;
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
  std::vector<std::any> splines;
  for (int i = 0; i < directoryEntrySection.size() * 0.5; ++i) {
    int entityType = std::stoi(trim(directoryEntrySection[i * 2].substr(5, 3)));
    if (entityType == 126) {
      splines.push_back(Create1DSpline(ParameterSectionToVector(parameterDataSection,
                                                                GetParameterSectionStartEndPointers(
                                                                    directoryEntrySection, i))));
    } else if (entityType == 128) {
      splines.push_back(Create2DSpline(ParameterSectionToVector(parameterDataSection,
                                                                GetParameterSectionStartEndPointers(
                                                                    directoryEntrySection, i))));
    }
  }
  return splines;
}

std::any io::IGESReader::Create1DSpline(const std::vector<double> &parameterData) {
  int upperSumIndex = static_cast<int>(parameterData[1]);
  std::array<Degree, 1> degree{};
  degree[0] = Degree{static_cast<int>(parameterData[2])};
  KnotVectors<1> knot_vector;
  std::vector<double> weights;
  std::vector<baf::ControlPoint> control_points;
  std::array<int, 2> knotsStartEnd{};
  std::array<int, 2> weightsStartEnd{};
  std::array<int, 2> controlPointsStartEnd{};
  knotsStartEnd[0] = 7;
  knotsStartEnd[1] = knotsStartEnd[0] + upperSumIndex + degree[0].get() + 1;
  weightsStartEnd[0] = knotsStartEnd[1] + 1;
  weightsStartEnd[1] = weightsStartEnd[0] + upperSumIndex;
  controlPointsStartEnd[0] = weightsStartEnd[1] + 1;
  controlPointsStartEnd[1] = controlPointsStartEnd[0] + (3 * upperSumIndex) + 2;
  std::vector<ParamCoord> knots;
  for (int i = knotsStartEnd[0]; i <= knotsStartEnd[1]; ++i) {
    knots.emplace_back(parameterData[i]);
  }
  knot_vector[0] = std::make_shared<baf::KnotVector>(knots);
  for (int i = weightsStartEnd[0]; i <= weightsStartEnd[1]; ++i) {
    weights.push_back(parameterData[i]);
  }
  std::vector<double> controlPointCoordinates;
  for (int i = controlPointsStartEnd[0]; i <= controlPointsStartEnd[1]; ++i) {
    controlPointCoordinates.push_back(parameterData[i]);
  }
  for (auto i = 0u; i < controlPointCoordinates.size(); i += 3) {
    control_points.push_back(baf::ControlPoint({controlPointCoordinates[i],
                                                controlPointCoordinates[i + 1],
                                                controlPointCoordinates[i + 2]}));
  }
  std::array<int, 1> number_of_points{};
  for (int i = 0; i < 1; ++i) {
    number_of_points[i] = static_cast<int>(knot_vector[i]->GetNumberOfKnots() - degree[i].get() - 1);
  }
  if (parameterData[5] == 1) {
    auto spl = std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points);
    return std::make_any<std::shared_ptr<spl::BSpline<1>>>(spl);
  } else {
    auto spl = std::make_shared<spl::NURBS<1>>(knot_vector, degree, control_points, weights);
    return std::make_any<std::shared_ptr<spl::NURBS<1>>>(spl);
  }
}

std::any io::IGESReader::Create2DSpline(const std::vector<double> &parameterData) {
  std::array<int, 2> upperSumIndex{};
  upperSumIndex[0] = static_cast<int>(parameterData[1]);
  upperSumIndex[1] = static_cast<int>(parameterData[2]);
  std::array<Degree, 2> degree{};
  degree[0] = Degree{static_cast<int>(parameterData[3])};
  degree[1] = Degree{static_cast<int>(parameterData[4])};
  KnotVectors<2> knot_vector;
  std::vector<double> weights;
  std::vector<baf::ControlPoint> control_points;
  std::array<std::array<int, 2>, 2> knotsStartEnd{};
  std::array<int, 2> weightsStartEnd{};
  std::array<int, 2> controlPointsStartEnd{};
  knotsStartEnd[0][0] = 10;
  knotsStartEnd[0][1] = knotsStartEnd[0][0] + upperSumIndex[0] + degree[0].get() + 1;
  knotsStartEnd[1][0] = knotsStartEnd[0][1] + 1;
  knotsStartEnd[1][1] = knotsStartEnd[1][0] + upperSumIndex[1] + degree[1].get() + 1;
  weightsStartEnd[0] = knotsStartEnd[1][1] + 1;
  weightsStartEnd[1] = weightsStartEnd[0] - 1 + ((1 + upperSumIndex[0]) * (1 + upperSumIndex[1]));
  controlPointsStartEnd[0] = weightsStartEnd[1] + 1;
  controlPointsStartEnd[1] = controlPointsStartEnd[0] - 1 + (3 * (1 + upperSumIndex[0]) * (1 + upperSumIndex[1])) + 2;
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
  for (auto i = 0u; i < controlPointCoordinates.size() - 2; i += 3) {
    control_points.push_back(baf::ControlPoint({controlPointCoordinates[i],
                                                controlPointCoordinates[i + 1],
                                                controlPointCoordinates[i + 2]}));
  }
  std::array<int, 2> number_of_points{};
  for (int i = 0; i < 2; ++i) {
    number_of_points[i] = static_cast<int>(knot_vector[i]->GetNumberOfKnots() - degree[i].get() - 1);
  }
  if (parameterData[7] == 1) {
    auto spl = std::make_shared<spl::BSpline<2>>(knot_vector, degree, control_points);
    return std::make_any<std::shared_ptr<spl::BSpline<2>>>(spl);
  } else {
    auto spl = std::make_shared<spl::NURBS<2>>(knot_vector, degree, control_points, weights);
    return std::make_any<std::shared_ptr<spl::NURBS<2>>>(spl);
  }
}

std::array<int, 2> io::IGESReader::GetParameterSectionStartEndPointers(std::vector<std::string> directoryEntrySection,
                                                                       int entityToBeRead) {
  std::string parameterDataStartPointer = trim(directoryEntrySection[entityToBeRead * 2].substr(8, 8));
  std::string parameterDataLineCount = trim(directoryEntrySection[entityToBeRead * 2 + 1].substr(24, 8));
  std::array<int, 2> ParameterSectionStartEndPointers{};
  ParameterSectionStartEndPointers[0] = std::stoi(trim(parameterDataStartPointer));
  ParameterSectionStartEndPointers[1] =
      std::stoi(trim(parameterDataStartPointer)) + std::stoi(trim(parameterDataLineCount)) - 1;
  return ParameterSectionStartEndPointers;
}

std::vector<double> io::IGESReader::ParameterSectionToVector(std::vector<std::string> parameterSection,
                                                             std::array<int, 2> ParameterSectionStartEndPointers) {
  int first = ParameterSectionStartEndPointers[0] - 1;
  int last = ParameterSectionStartEndPointers[1] - 1;
  std::string temp;
  for (int i = first; i <= last; ++i) {
    temp.append(parameterSection[i]);
  }
  return util::StringOperations::DelimitedStringToVector(temp);
}

std::string io::IGESReader::trim(std::string s) {
  return util::StringOperations::trim(std::move(s));
}
