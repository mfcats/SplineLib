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

namespace io {
template<int DIM>
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
      vector_of_splines.emplace_back(GetSpline(spline_position, entries));
    }
    return vector_of_splines;
  }

 private:
  std::vector<int> GetSplinePositions(const std::vector<std::string> &entries) const {
    std::vector<int> spline_positions;
    for (int i = 0; i < entries.size(); i++) {
      if (entries[i] == "BSPLINE" && FitsDimension(i, entries)) spline_positions.push_back(i - 1);
    }
    return spline_positions;
  }

  std::any GetSpline(int start, const std::vector<std::string> &entries) const {
    std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vector = GetKnotVectors(start, entries);
    std::array<Degree, DIM> degree = GetDegrees(start, entries);
    bool rational = IsRational(start, entries);
    std::vector<baf::ControlPoint> control_points = GetControlPoints(start, entries, rational);
    if (!rational) {
      return std::make_any<std::shared_ptr<spl::BSpline<DIM>>>(std::make_shared<spl::BSpline<DIM>>(knot_vector,
                                                                                                   degree,
                                                                                                   control_points));
    } else {
      std::vector<double> weights = GetWeights(start, entries);
      return std::make_any<std::shared_ptr<spl::NURBS<DIM>>>(std::make_shared<spl::NURBS<DIM>>(knot_vector,
                                                                                               degree,
                                                                                               control_points,
                                                                                               weights));
    }
  }

  std::array<Degree, DIM> GetDegrees(int start, const std::vector<std::string> &entries) const {
    std::array<Degree, DIM> degrees;
    for (int i = 0; i < DIM; i++) {
      degrees[i] =
          Degree(util::StringOperations::StringVectorToNumberVector<int>({entries[start + DIM + 2 + i]})[0] - 1);
    }
    return degrees;
  }

  std::array<std::shared_ptr<baf::KnotVector>, DIM>
  GetKnotVectors(int start, const std::vector<std::string> &entries) const {
    std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vectors;
    for (int i = 0; i < DIM; i++) {
      while (!StartsWith(entries[start++], "[KV")) {}
      std::vector<ParamCoord> knots;
      while (!StartsWith(entries[start], "[")) {
        knots.emplace_back(util::StringOperations::StringVectorToNumberVector<double>({entries[start++]})[0]);
      }
      knot_vectors[i] = std::make_shared<baf::KnotVector>(knots);
    }
    return knot_vectors;
  }

  std::vector<baf::ControlPoint>
  GetControlPoints(int start, const std::vector<std::string> &entries, bool rational) const {
    std::vector<baf::ControlPoint> control_points;
    int number_of_control_points = GetNumberOfControlPoints(start, entries);
    start = GetPositionOfFirstControlPoint(start, entries);
    for (int i = 0; i < number_of_control_points; i++) {
      std::vector<double> coordinates;
      while (!EndsWith(entries[start], "]")) {
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
      while (!EndsWith(entries[start++], "]")) {}
    }
    return weights;
  }

  int GetNumberOfControlPoints(int start, const std::vector<std::string> &entries) const {
    int total_number_of_points = 1;
    for (int i = 0; i < DIM; i++) {
      total_number_of_points *= util::StringOperations::StringVectorToNumberVector<int>({entries[start + 2 + i]})[0];
    }
    return total_number_of_points;
  }

  int GetPositionOfFirstControlPoint(int start, const std::vector<std::string> &entries) const {
    start++;
    while (!StartsWith(entries[start], "[") || StartsWith(entries[start], "[KV")) {
      start++;
    }
    return start;
  }

  bool StartsWith(const std::string &string, const std::string &start_of_string) const {
    return string.find(start_of_string) == 0;
  }

  bool EndsWith(const std::string &string, const std::string &end_of_string) const {
    return string.find(end_of_string) == string.length() - end_of_string.length();
  }

  std::string trim(const std::string &string) const {
    if (StartsWith(string, "[")) {
      return string.substr(1, string.length() - 1);
    } else if (EndsWith(string, "]")) {
      return string.substr(0, string.length() - 1);
    } else {
      return string;
    }
  }

  bool FitsDimension(int position, const std::vector<std::string> &entries) const {
    return ((DIM == 1 && entries[position - 1] == "[CURVE") || (DIM == 2 && entries[position - 1] == "[SURFACE")
        || (DIM == 3 && entries[position - 1] == "[TRIVAR"));
  }

  bool IsRational(int start_of_spline, const std::vector<std::string> &entries) const {
    return StartsWith(entries[start_of_spline + 2 * DIM + 2], "P");
  }
};
}  // namespace io

#endif  // SRC_IO_IRIT_READER_H_
