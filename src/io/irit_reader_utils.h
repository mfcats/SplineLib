/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_IRIT_READER_UTILS_H_
#define SRC_IO_IRIT_READER_UTILS_H_

#include <string>
#include <vector>

#include "b_spline.h"
#include "nurbs.h"
#include "string_operations.h"

namespace io {
class IRITReaderUtils {
 public:
  template<int DIM>
  static std::array<Degree, DIM> GetDegrees(int start, const std::vector<std::string> &entries) {
    std::array<Degree, DIM> degrees;
    for (int i = 0; i < DIM; i++) {
      degrees[i] =
          Degree(util::StringOperations::StringVectorToNumberVector<int>({entries[start + DIM + 2 + i]})[0] - 1);
    }
    return degrees;
  }

  template<int DIM>
  static KnotVectors<DIM> GetKnotVectors(int start, const std::vector<std::string> &entries) {
    KnotVectors<DIM> knot_vectors;
    for (int i = 0; i < DIM; i++) {
      while (!util::StringOperations::StartsWith(entries[start++], "[KV")) {}
      std::vector<ParamCoord> knots;
      while (!util::StringOperations::StartsWith(entries[start], "[")) {
        knots.emplace_back(util::StringOperations::StringToDouble(util::StringOperations::trim(entries[start++])));
      }
      knot_vectors[i] = std::make_shared<baf::KnotVector>(knots);
    }
    return knot_vectors;
  }

  template<int DIM>
  static bool IsRational(int start_of_spline, const std::vector<std::string> &entries) {
    return util::StringOperations::StartsWith(entries[start_of_spline + 2 * DIM + 2], "P");
  }
};
}  // namespace io

#endif  // SRC_IO_IRIT_READER_UTILS_H_
