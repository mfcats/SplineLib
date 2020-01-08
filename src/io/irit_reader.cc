/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "src/io/irit_reader.h"

#include <fstream>

#include "src/io/irit_reader_utils.h"
#include "src/spl/nurbs.h"
#include "src/util/string_operations.h"

namespace splinelib::src::io {
std::vector<std::any> IRITReader::ReadFile(const char *filename) {
  std::vector<std::any> vector_of_splines;
  std::ifstream newFile;
  newFile.open(filename);
  if (!newFile.good()) {
    throw std::runtime_error("IRIT file could not be opened.");
  }
  std::string line;
  std::string file;
  while (getline(newFile, line)) {
    file += line;
  }
  std::vector<std::string> entries = util::string_operations::SplitStringAtDelimiter(file, ' ');
  std::vector<int> spline_positions = GetSplinePositions(entries);
  for (int &spline_position : spline_positions) {
    if (GetDimension(entries[spline_position]) == 1) {
      vector_of_splines.emplace_back(Get1DSpline(spline_position, entries));
    } else if (GetDimension(entries[spline_position]) == 2) {
      vector_of_splines.emplace_back(Get2DSpline(spline_position, entries));
    } else if (GetDimension(entries[spline_position]) == 3) {
      vector_of_splines.emplace_back(Get3DSpline(spline_position, entries));
    }
  }
  return vector_of_splines;
}

std::vector<int> IRITReader::GetSplinePositions(const std::vector<std::string> &entries) const {
  std::vector<int> spline_positions;
  for (auto i = 0u; i < entries.size(); i++) {
    if (entries[i] == "BSPLINE") {
      spline_positions.push_back(i - 1);
    }
  }
  return spline_positions;
}

int IRITReader::GetDimension(const std::string &type) {
  return type == "[CURVE" ? 1 : (type == "[SURFACE" ? 2 : (type == "[TRIVAR" ? 3 : 0));
}

std::any IRITReader::Get1DSpline(int start, const std::vector<std::string> &entries) const {
  baf::KnotVectors<1> knot_vector = IRITReaderUtils::GetKnotVectors<1>(start, entries);
  std::array<Degree, 1> degree = IRITReaderUtils::GetDegrees<1>(start, entries);
  bool rational = IRITReaderUtils::IsRational<1>(start, entries);
  auto weights = GetWeights(start, entries, rational);
  std::vector<spl::ControlPoint> control_points = GetControlPoints(start, entries, rational, weights);
  if (!rational) {
    return std::make_any<std::shared_ptr<spl::BSpline<1>>>(
        std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points));
  }
  return std::make_any<std::shared_ptr<spl::NURBS<1>>>(
      std::make_shared<spl::NURBS<1>>(knot_vector, degree, control_points, weights));
}

std::any IRITReader::Get2DSpline(int start, const std::vector<std::string> &entries) const {
  baf::KnotVectors<2> knot_vector = IRITReaderUtils::GetKnotVectors<2>(start, entries);
  std::array<Degree, 2> degree = IRITReaderUtils::GetDegrees<2>(start, entries);
  bool rational = IRITReaderUtils::IsRational<2>(start, entries);
  auto weights = GetWeights(start, entries, rational);
  std::vector<spl::ControlPoint> control_points = GetControlPoints(start, entries, rational, weights);
  if (!rational) {
    return std::make_any<std::shared_ptr<spl::BSpline<2>>>(
        std::make_shared<spl::BSpline<2>>(knot_vector, degree, control_points));
  }
  return std::make_any<std::shared_ptr<spl::NURBS<2>>>(
      std::make_shared<spl::NURBS<2>>(knot_vector, degree, control_points, weights));
}

std::any IRITReader::Get3DSpline(int start, const std::vector<std::string> &entries) const {
  baf::KnotVectors<3> knot_vector = IRITReaderUtils::GetKnotVectors<3>(start, entries);
  std::array<Degree, 3> degree = IRITReaderUtils::GetDegrees<3>(start, entries);
  bool rational = IRITReaderUtils::IsRational<3>(start, entries);
  auto weights = GetWeights(start, entries, rational);
  std::vector<spl::ControlPoint> control_points = GetControlPoints(start, entries, rational, weights);
  if (!rational) {
    return std::make_any<std::shared_ptr<spl::BSpline<3>>>(
        std::make_shared<spl::BSpline<3>>(knot_vector, degree, control_points));
  }
  return std::make_any<std::shared_ptr<spl::NURBS<3>>>(
      std::make_shared<spl::NURBS<3>>(knot_vector, degree, control_points, weights));
}

int IRITReader::GetTotalNumberOfControlPoints(int start, const std::vector<std::string> &entries) {
  int total_number_of_points = 1;
  for (int i = 0; i < GetDimension(entries[start]); i++) {
    total_number_of_points *=
        util::string_operations::ConvertStringVectorToNumberVector<int>({entries[start + 2 + i]})[0];
  }
  return total_number_of_points;
}

std::vector<spl::ControlPoint> IRITReader::GetControlPoints(int start,
                                                                const std::vector<std::string> &entries,
                                                                bool rational,
                                                                std::vector<double> weights) const {
  std::vector<spl::ControlPoint> control_points;
  int number_of_control_points = GetTotalNumberOfControlPoints(start, entries);
  start = GetPositionOfFirstControlPoint(start, entries);
  for (int i = 0; i < number_of_control_points; i++) {
    std::vector<double> coordinates;
    while (!util::string_operations::EndsWith(entries[start], "]")) {
      coordinates.push_back(
          util::string_operations::ConvertStringToNumber<double>(util::string_operations::TrimSpacesAndSquareBrackets(
              entries[start++])) /
          weights[i]);
    }
    coordinates.push_back(
        util::string_operations::ConvertStringToNumber<double>(util::string_operations::TrimSpacesAndSquareBrackets(
            entries[start++])) / weights[i]);
    if (rational) {
      coordinates.erase(coordinates.begin());
    }
    control_points.emplace_back(coordinates);
  }
  return control_points;
}

std::vector<double> IRITReader::GetWeights(int start,
                                               const std::vector<std::string> &entries,
                                               bool rational) const {
  auto number_of_control_points = static_cast<size_t>(GetTotalNumberOfControlPoints(start, entries));
  std::vector<double> weights(number_of_control_points, 1.0);
  if (rational) {
    start = GetPositionOfFirstControlPoint(start, entries);
    for (auto i = 0u; i < number_of_control_points; i++) {
      weights[i] =
          util::string_operations::ConvertStringToNumber<double>(util::string_operations::TrimSpacesAndSquareBrackets(
              entries[start++]));
      while (!util::string_operations::EndsWith(entries[start++], "]")) {}
    }
  }
  return weights;
}

int IRITReader::GetPositionOfFirstControlPoint(int start, const std::vector<std::string> &entries) const {
  ++start;
  while (!util::string_operations::StartsWith(entries[start], "[")
      || util::string_operations::StartsWith(entries[start], "[KV")) {
    ++start;
  }
  return start;
}
}  // namespace splinelib::src::io
