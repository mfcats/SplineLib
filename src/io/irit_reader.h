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
    for (int i = 0; i < spline_positions.size(); i += 2) {
      vector_of_splines.emplace_back(GetSpline(spline_positions[i], spline_positions[i + 1], entries));
    }
    return vector_of_splines;
  }

 private:
  std::vector<int> GetSplinePositions(const std::vector<std::string> &entries) const {
    std::vector<int> spline_positions;
    for (int i = 0; i < entries.size(); i++) {
      if (entries[i] == "BSPLINE" && FitsDimension(i, entries)) {
        spline_positions.push_back(i - 1);
        int position = i - 1;
        int brace_difference = 0;
        for (; position < entries.size(); position++) {
          if (StartsWith(entries[position], "[")) brace_difference++;
          else if (EndsWith(entries[position], "]")) brace_difference--;
          if (brace_difference == 0) break;
        }
        spline_positions.push_back(position);
      }
    }
    return spline_positions;
  }

  std::any GetSpline(int start_of_spline, int end_of_spline, std::vector<std::string> entries) const {
    std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vector =
        GetKnotVectors(start_of_spline, end_of_spline, entries);
    std::array<Degree, DIM> degree = GetDegrees(start_of_spline, entries);
    std::vector<baf::ControlPoint> control_points =
        GetControlPoints(start_of_spline, entries);
    return std::make_any<spl::BSpline<DIM>>(knot_vector, degree, control_points);
  }

  std::array<Degree, DIM> GetDegrees(int start_of_spline, const std::vector<std::string> &entries) const {
    std::array<Degree, DIM> degrees;
    for (int i = 0; i < DIM; i++) {
      degrees[i] = Degree(
          util::StringOperations::StringVectorToNumberVector<int>({entries[start_of_spline + DIM + 2 + i]})[0] - 1);
    }
    return degrees;
  }

  std::array<std::shared_ptr<baf::KnotVector>, DIM>
  GetKnotVectors(int start_of_spline, int end_of_spline, const std::vector<std::string> &entries) const {
    std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vectors;
    int count = 0;
    for (int i = start_of_spline; i < end_of_spline; i++) {
      if (StartsWith(entries[i], "[KV")) {
        std::vector<ParamCoord> knots;
        i++;
        while (!StartsWith(entries[i], "[")) {
          knots.emplace_back(util::StringOperations::StringVectorToNumberVector<double>({entries[i]})[0]);
          i++;
        }
        knot_vectors[count] = std::make_shared<baf::KnotVector>(knots);
        count++;
        i--;
      }
    }
    return knot_vectors;
  }

  std::vector<baf::ControlPoint>
  GetControlPoints(int start_of_spline, const std::vector<std::string> &entries) const {
    std::vector<baf::ControlPoint> control_points;
    std::array<int, DIM> number_of_points;
    int total_number_of_points = 1;
    for (int i = 0; i < DIM; i++) {
      number_of_points[i] =
          util::StringOperations::StringVectorToNumberVector<int>({entries[start_of_spline + 2 + i]})[0];
      total_number_of_points *= number_of_points[i];
    }
    start_of_spline++;
    while (!StartsWith(entries[start_of_spline], "[") || StartsWith(entries[start_of_spline], "[KV")) {
      start_of_spline++;
    }
    for (int i = 0; i < total_number_of_points; i++) {
      std::vector<double> coordinates;
      while (!EndsWith(entries[start_of_spline], "]")) {
        coordinates.push_back(
            util::StringOperations::StringVectorToNumberVector<double>({trim(entries[start_of_spline])})[0]);
        start_of_spline++;
      }
      coordinates.push_back(util::StringOperations::StringVectorToNumberVector<double>(
          {trim(entries[start_of_spline])})[0]);
      control_points.emplace_back(coordinates);
      start_of_spline++;
    }
    return control_points;
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
    return ((DIM == 1 && entries[position - 1] == "[CURVE") || (DIM == 2 && entries[position - 1] == "[SURFACE"));
  }
};
}  // namespace io

#endif  // SRC_IO_IRIT_READER_H_
