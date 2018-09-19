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
        spline_positions.push_back(ClosingBrace(i - 1, entries));
      }
    }
    return spline_positions;
  }

  std::any GetSpline(int first, int last, const std::vector<std::string> &entries) const {
    std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vector = GetKnotVectors(first, last, entries);
    std::array<Degree, DIM> degree = GetDegrees(first, entries);
    bool rational = IsRational(first, entries);
    std::vector<baf::ControlPoint> control_points = GetControlPoints(first, entries, rational);
    if (!rational) {
      return std::make_any<spl::BSpline<DIM>>(knot_vector, degree, control_points);
    } else {
      std::vector<double> weights = GetWeights(first, entries);
      return std::make_any<spl::NURBS<DIM>>(knot_vector, degree, control_points, weights);
    }
  }

  std::array<Degree, DIM> GetDegrees(int first, const std::vector<std::string> &entries) const {
    std::array<Degree, DIM> degrees;
    for (int i = 0; i < DIM; i++) {
      degrees[i] =
          Degree(util::StringOperations::StringVectorToNumberVector<int>({entries[first + DIM + 2 + i]})[0] - 1);
    }
    return degrees;
  }

  std::array<std::shared_ptr<baf::KnotVector>, DIM>
  GetKnotVectors(int first, int last, const std::vector<std::string> &entries) const {
    std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vectors;
    for (int i = 0; i < DIM; i++) {
      while (!StartsWith(entries[first++], "[KV")) {}
      std::vector<ParamCoord> knots;
      while (!StartsWith(entries[first], "[")) {
        knots.emplace_back(util::StringOperations::StringVectorToNumberVector<double>({entries[first++]})[0]);
      }
      knot_vectors[i] = std::make_shared<baf::KnotVector>(knots);
    }
    return knot_vectors;
  }

  std::vector<baf::ControlPoint>
  GetControlPoints(int first, const std::vector<std::string> &entries, bool rational) const {
    std::vector<baf::ControlPoint> control_points;
    int number_of_control_points = GetNumberOfControlPoints(first, entries);
    first = GetPositionOfFirstControlPoint(first, entries);
    for (int i = 0; i < number_of_control_points; i++) {
      std::vector<double> coordinates;
      while (!EndsWith(entries[first], "]")) {
        coordinates.push_back(util::StringOperations::StringVectorToNumberVector<double>({trim(entries[first++])})[0]);
      }
      coordinates.push_back(util::StringOperations::StringVectorToNumberVector<double>({trim(entries[first++])})[0]);
      if (rational) coordinates.erase(coordinates.begin());
      control_points.emplace_back(coordinates);
    }
    return control_points;
  }

  std::vector<double> GetWeights(int first, const std::vector<std::string> &entries) const {
    std::vector<double> weights;
    int number_of_control_points = GetNumberOfControlPoints(first, entries);
    first = GetPositionOfFirstControlPoint(first, entries);
    for (int i = 0; i < number_of_control_points; i++) {
      weights.push_back(util::StringOperations::StringVectorToNumberVector<double>({trim(entries[first++])})[0]);
      while (!EndsWith(entries[first], "]")) {
        first++;
      }
      first++;
    }
    return weights;
  }

  int GetNumberOfControlPoints(int first, const std::vector<std::string> &entries) const {
    int total_number_of_points = 1;
    for (int i = 0; i < DIM; i++) {
      total_number_of_points *= util::StringOperations::StringVectorToNumberVector<int>({entries[first + 2 + i]})[0];
    }
    return total_number_of_points;
  }

  int GetPositionOfFirstControlPoint(int first, const std::vector<std::string> &entries) const {
    first++;
    while (!StartsWith(entries[first], "[") || StartsWith(entries[first], "[KV")) {
      first++;
    }
    return first;
  }

  int ClosingBrace(int position, const std::vector<std::string> &entries) const {
    int brace_difference = 0;
    for (; position < entries.size(); position++) {
      if (StartsWith(entries[position], "[")) brace_difference++;
      else if (EndsWith(entries[position], "]")) brace_difference--;
      if (brace_difference == 0) break;
    }
    return position;
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

  bool IsRational(int start_of_spline, const std::vector<std::string> &entries) const {
    return StartsWith(entries[start_of_spline + 2 * DIM + 2], "P");
  }
};
}  // namespace io

#endif  // SRC_IO_IRIT_READER_H_
