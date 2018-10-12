/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_IRIT_READER_H_
#define SRC_IO_IRIT_READER_H_

#include <any>
#include <fstream>
#include <string>
#include <vector>

#include "b_spline.h"
#include "nurbs.h"
#include "string_operations.h"
#include "irit_reader_utils.h"

namespace io {
class IRITReader {
 public:
  IRITReader() = default;

  std::vector<std::any> ReadIRITFile(const char *filename) {
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
    std::vector<std::string> entries = util::StringOperations::split(file, ' ');
    std::vector<int> spline_positions = GetSplinePositions(entries);
    for (int &spline_position : spline_positions) {
      switch (GetDimension(entries[spline_position])) {
        case 1: {
          vector_of_splines.emplace_back(Get1DSpline(spline_position, entries));
          break;
        }
        case 2: {
          vector_of_splines.emplace_back(Get2DSpline(spline_position, entries));
          break;
        }
        case 3: {
          vector_of_splines.emplace_back(Get3DSpline(spline_position, entries));
          break;
        }
        default: {}
      }
    }
    return vector_of_splines;
  }

 private:
  std::vector<int> GetSplinePositions(const std::vector<std::string> &entries) const {
    std::vector<int> spline_positions;
    for (auto i = 0u; i < entries.size(); i++) {
      if (entries[i] == "BSPLINE") spline_positions.push_back(i - 1);
    }
    return spline_positions;
  }

  static int GetDimension(const std::string &type) {
    return type == "[CURVE" ? 1 : (type == "[SURFACE" ? 2 : (type == "[TRIVAR" ? 3 : 0));
  }

  std::any Get1DSpline(int start, const std::vector<std::string> &entries) const {
    std::array<std::shared_ptr<baf::KnotVector>, 1>
        knot_vector = io::IRITReaderUtils<1>::GetKnotVectors(start, entries);
    std::array<Degree, 1> degree = io::IRITReaderUtils<1>::GetDegrees(start, entries);
    bool rational = io::IRITReaderUtils<1>::IsRational(start, entries);
    std::vector<baf::ControlPoint> control_points = GetControlPoints(start, entries, rational);
    if (!rational) {
      return std::make_any<std::shared_ptr<spl::BSpline<1>>>(
          std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points));
    } else {
      return std::make_any<std::shared_ptr<spl::NURBS<1>>>(
          std::make_shared<spl::NURBS<1>>(knot_vector, degree, control_points, GetWeights(start, entries)));
    }
  }

  std::any Get2DSpline(int start, const std::vector<std::string> &entries) const {
    std::array<std::shared_ptr<baf::KnotVector>, 2>
        knot_vector = io::IRITReaderUtils<2>::GetKnotVectors(start, entries);
    std::array<Degree, 2> degree = io::IRITReaderUtils<2>::GetDegrees(start, entries);
    bool rational = io::IRITReaderUtils<2>::IsRational(start, entries);
    std::vector<baf::ControlPoint> control_points = GetControlPoints(start, entries, rational);
    if (!rational) {
      return std::make_any<std::shared_ptr<spl::BSpline<2>>>(
          std::make_shared<spl::BSpline<2>>(knot_vector, degree, control_points));
    } else {
      return std::make_any<std::shared_ptr<spl::NURBS<2>>>(
          std::make_shared<spl::NURBS<2>>(knot_vector, degree, control_points, GetWeights(start, entries)));
    }
  }

  std::any Get3DSpline(int start, const std::vector<std::string> &entries) const {
    std::array<std::shared_ptr<baf::KnotVector>, 3>
        knot_vector = io::IRITReaderUtils<3>::GetKnotVectors(start, entries);
    std::array<Degree, 3> degree = io::IRITReaderUtils<3>::GetDegrees(start, entries);
    bool rational = io::IRITReaderUtils<3>::IsRational(start, entries);
    std::vector<baf::ControlPoint> control_points = GetControlPoints(start, entries, rational);
    if (!rational) {
      return std::make_any<std::shared_ptr<spl::BSpline<3>>>(
          std::make_shared<spl::BSpline<3>>(knot_vector, degree, control_points));
    } else {
      return std::make_any<std::shared_ptr<spl::NURBS<3>>>(
          std::make_shared<spl::NURBS<3>>(knot_vector, degree, control_points, GetWeights(start, entries)));
    }
  }

  static int GetNumberOfControlPoints(int start, const std::vector<std::string> &entries) {
    int total_number_of_points = 1;
    for (int i = 0; i < GetDimension(entries[start]); i++) {
      total_number_of_points *= util::StringOperations::StringVectorToNumberVector<int>({entries[start + 2 + i]})[0];
    }
    return total_number_of_points;
  }

  std::vector<baf::ControlPoint>
  GetControlPoints(int start, const std::vector<std::string> &entries, bool rational) const {
    std::vector<baf::ControlPoint> control_points;
    int number_of_control_points = GetNumberOfControlPoints(start, entries);
    start = GetPositionOfFirstControlPoint(start, entries);
    for (int i = 0; i < number_of_control_points; i++) {
      std::vector<double> coordinates;
      while (!util::StringOperations::EndsWith(entries[start], "]")) {
        coordinates.push_back(util::StringOperations::StringVectorToNumberVector<double>({trim(entries[start++])})[0]);
      }
      coordinates.push_back(util::StringOperations::StringVectorToNumberVector<double>({trim(entries[start++])})[0]);
      if (rational) coordinates.erase(coordinates.begin());
      control_points.emplace_back(coordinates);
    }
    return control_points;
  }

  std::vector<double> GetWeights(int start, const std::vector<std::string> &entries) const {
    std::vector<double> weights;
    int number_of_control_points = GetNumberOfControlPoints(start, entries);
    start = GetPositionOfFirstControlPoint(start, entries);
    for (int i = 0; i < number_of_control_points; i++) {
      weights.push_back(util::StringOperations::StringVectorToNumberVector<double>({trim(entries[start++])})[0]);
      while (!util::StringOperations::EndsWith(entries[start++], "]")) {}
    }
    return weights;
  }

  int GetPositionOfFirstControlPoint(int start, const std::vector<std::string> &entries) const {
    start++;
    while (!util::StringOperations::StartsWith(entries[start], "[")
        || util::StringOperations::StartsWith(entries[start], "[KV")) {
      start++;
    }
    return start;
  }

  std::string trim(const std::string &string) const {
    if (util::StringOperations::StartsWith(string, "[")) {
      return string.substr(1, string.length() - 1);
    } else if (util::StringOperations::EndsWith(string, "]")) {
      return string.substr(0, string.length() - 1);
    } else {
      return string;
    }
  }
};
}  // namespace io

#endif  // SRC_IO_IRIT_READER_H_
