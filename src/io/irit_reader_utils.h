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
#include "src/util/string_operations.h"

namespace splinelib::src::io {
class IRITReaderUtils {
 public:
  template<int PARAMETRIC_DIMENSIONALITY>
  static std::array<Degree, PARAMETRIC_DIMENSIONALITY> GetDegrees(int start, const std::vector<std::string> &entries) {
    std::array<Degree, PARAMETRIC_DIMENSIONALITY> degrees{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; i++) {
      degrees[i] =
          Degree(util::string_operations::StringVectorToNumberVector<int>({entries[start + PARAMETRIC_DIMENSIONALITY + 2
              + i]})[0] - 1);
    }
    return degrees;
  }

  template<int PARAMETRIC_DIMENSIONALITY>
  static baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> GetKnotVectors(int start,
                                                                    const std::vector<std::string> &entries) {
    baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> knot_vectors;
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; i++) {
      while (!util::string_operations::StartsWith(entries[start++], "[KV")) {}
      std::vector<ParametricCoordinate> knots;
      while (!util::string_operations::StartsWith(entries[start], "[")) {
        knots.emplace_back(
            util::string_operations::StringToNumber<double>(util::string_operations::Trim(entries[start++])));
      }
      knot_vectors[i] = std::make_shared<baf::KnotVector>(knots);
    }
    return knot_vectors;
  }

  template<int PARAMETRIC_DIMENSIONALITY>
  static bool IsRational(int start_of_spline, const std::vector<std::string> &entries) {
    return util::string_operations::StartsWith(entries[start_of_spline + 2 * PARAMETRIC_DIMENSIONALITY + 2], "P");
  }
};
}  // namespace splinelib::src::io

#endif  // SRC_IO_IRIT_READER_UTILS_H_
