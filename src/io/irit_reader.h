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

  std::vector<std::any> ReadFile(const char *filename) {
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
    KnotVectors<1> knot_vector = io::IRITReaderUtils::GetKnotVectors<1>(start, entries);
    std::array<Degree, 1> degree = io::IRITReaderUtils::GetDegrees<1>(start, entries);
    bool rational = io::IRITReaderUtils::IsRational<1>(start, entries);
    std::vector<baf::ControlPoint> control_points = GetControlPoints(start, entries, rational);
    if (!rational) {
      return std::make_any<std::shared_ptr<spl::BSpline<1>>>(
          std::make_shared<spl::BSpline<1>>(knot_vector, degree, control_points));
    } else {
      auto weights = GetWeights(start, entries);
      return std::make_any<std::shared_ptr<spl::NURBS<1>>>(
          std::make_shared<spl::NURBS<1>>(knot_vector, degree,
                                          GetEuclidianControlPoints(control_points, weights), weights));
    }
  }

  std::any Get2DSpline(int start, const std::vector<std::string> &entries) const {
    KnotVectors<2> knot_vector = io::IRITReaderUtils::GetKnotVectors<2>(start, entries);
    std::array<Degree, 2> degree = io::IRITReaderUtils::GetDegrees<2>(start, entries);
    bool rational = io::IRITReaderUtils::IsRational<2>(start, entries);
    std::vector<baf::ControlPoint> control_points = GetControlPoints(start, entries, rational);
    if (!rational) {
      return std::make_any<std::shared_ptr<spl::BSpline<2>>>(
          std::make_shared<spl::BSpline<2>>(knot_vector, degree, control_points));
    } else {
      auto weights = GetWeights(start, entries);
      return std::make_any<std::shared_ptr<spl::NURBS<2>>>(
          std::make_shared<spl::NURBS<2>>(knot_vector, degree,
                                          GetEuclidianControlPoints(control_points, weights), weights));
    }
  }

  std::any Get3DSpline(int start, const std::vector<std::string> &entries) const {
    KnotVectors<3> knot_vector = io::IRITReaderUtils::GetKnotVectors<3>(start, entries);
    std::array<Degree, 3> degree = io::IRITReaderUtils::GetDegrees<3>(start, entries);
    bool rational = io::IRITReaderUtils::IsRational<3>(start, entries);
    std::vector<baf::ControlPoint> control_points = GetControlPoints(start, entries, rational);
    if (!rational) {
      return std::make_any<std::shared_ptr<spl::BSpline<3>>>(
          std::make_shared<spl::BSpline<3>>(knot_vector, degree, control_points));
    } else {
      auto weights = GetWeights(start, entries);
      return std::make_any<std::shared_ptr<spl::NURBS<3>>>(
          std::make_shared<spl::NURBS<3>>(knot_vector, degree,
                                          GetEuclidianControlPoints(control_points, weights), weights));
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
        coordinates.push_back(util::StringOperations::StringToDouble(util::StringOperations::trim(entries[start++])));
      }
      coordinates.push_back(util::StringOperations::StringToDouble(util::StringOperations::trim(entries[start++])));
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
      weights.push_back(util::StringOperations::StringToDouble(util::StringOperations::trim(entries[start++])));
      while (!util::StringOperations::EndsWith(entries[start++], "]")) {}
    }
    return weights;
  }

  std::vector<baf::ControlPoint>
  GetEuclidianControlPoints(std::vector<baf::ControlPoint> homogeneous_points, std::vector<double> weights) const {
    std::vector<baf::ControlPoint> control_points;
    for (auto i = 0u; i < homogeneous_points.size(); ++i) {
      std::vector<double> coordinates;
      for (int j = 0; j < homogeneous_points[0].GetDimension(); ++j) {
        coordinates.push_back(homogeneous_points[i].GetValue(j) / weights[i]);
      }
      baf::ControlPoint point(coordinates);
      control_points.push_back(point);
    }
    return control_points;
  }

  int GetPositionOfFirstControlPoint(int start, const std::vector<std::string> &entries) const {
    ++start;
    while (!util::StringOperations::StartsWith(entries[start], "[")
        || util::StringOperations::StartsWith(entries[start], "[KV")) {
      ++start;
    }
    return start;
  }
};
}  // namespace io

#endif  // SRC_IO_IRIT_READER_H_
