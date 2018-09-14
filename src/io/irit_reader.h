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
#include <iostream>
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
    std::cout << std::endl;
    for (int i = 0; i < entries.size(); i++) {
      std::cout << i << ": " << entries[i] << std::endl;
    }

    std::vector<int> start_of_splines;
    for (int i = 0; i < entries.size(); i++) {
      if (entries[i] == "BSPLINE") {
        start_of_splines.push_back(i - 1);
      }
    }
    std::vector<int> end_of_splines;
    for (const auto &spline : start_of_splines) {
      end_of_splines.push_back(GetEndOfSpline(spline, entries));
    }

    for (int i = 0; i < start_of_splines.size(); i++) {
      std::cout << "start: " << start_of_splines[i] << ", end: " << end_of_splines[i] << std::endl;
    }

    for (int i = 0; i < start_of_splines.size(); i++) {
      std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vector =
          GetKnotVectors(start_of_splines[i], end_of_splines[i], entries);
      std::array<Degree, DIM> degree = GetDegrees(start_of_splines[i], end_of_splines[i], entries);
      std::vector<baf::ControlPoint> control_points =
          GetControlPoints(start_of_splines[i], end_of_splines[i], entries);
      for (auto cp : control_points) {
        for (int j = 0; j < (cp.GetDimension() - 1); j++) {
          std::cout << cp.GetValue(j) << ", ";
        }
        std::cout << cp.GetValue(cp.GetDimension() - 1) << std::endl;
      }
      vector_of_splines.emplace_back(std::make_any<spl::BSpline<DIM>>(knot_vector, degree, control_points));
    }
    return vector_of_splines;
  }

 private:
  int GetEndOfSpline(int start, const std::vector<std::string> &entries) const {
    int count = 0;
    int i = start;
    for (; i < entries.size(); i++) {
      if (StartsWith(entries[i], "[")) {
        count++;
      } else if (EndsWith(entries[i], "]")) {
        count--;
      }
      if (count == 0) {
        break;
      }
    }
    return i;
  }

  std::array<Degree, DIM>
  GetDegrees(int start_of_spline, int end_of_spline, const std::vector<std::string> &entries) {
    std::array<Degree, DIM> degrees;
    for (int i = 0; i < DIM; i++) {
      degrees[i] =
          Degree(
              util::StringOperations::StringVectorToNumberVector<int>({entries[start_of_spline + DIM + 2 + i]})[0] - 1);
    }
    return degrees;
  }

  std::array<std::shared_ptr<baf::KnotVector>, DIM>
  GetKnotVectors(int start_of_spline, int end_of_spline, const std::vector<std::string> &entries) {
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
      }
    }
    return knot_vectors;
  }

  std::vector<baf::ControlPoint>
  GetControlPoints(int start_of_spline, int end_of_spline, const std::vector<std::string> &entries) {
    std::vector<baf::ControlPoint> control_points;
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; i++) {
      number_of_points[i] =
          util::StringOperations::StringVectorToNumberVector<int>({entries[start_of_spline + 2 + i]})[0];
    }
    while (!StartsWith(entries[start_of_spline], "[KV")) {
      start_of_spline++;
    }
    while (!EndsWith(entries[start_of_spline], "]")) {
      start_of_spline++;
    }
    start_of_spline++;
    for (int i = 0; i < number_of_points[0]; i++) {
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
};
}  // namespace io

#endif  // SRC_IO_IRIT_READER_H_
