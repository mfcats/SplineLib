/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_XML_READER_UTILS_H_
#define SRC_IO_XML_READER_UTILS_H_

#include <string>
#include <vector>

#include "external/pugixml/pugixml.hpp"

#include "src/spl/b_spline.h"
#include "src/spl/nurbs.h"
#include "src/util/string_operations.h"

namespace splinelib::src::io {
template<int PARAMETRIC_DIMENSIONALITY>
class XMLReaderUtils {
 public:
  static std::array<Degree, PARAMETRIC_DIMENSIONALITY> GetDegrees(pugi::xml_node *spline) {
    return StringVectorToDegreeArray(
        util::string_operations::SplitStringAtDelimiter(spline->child("deg").first_child().value(), ' '));
  }

  static baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> GetKnotVectors(pugi::xml_node *spline) {
    baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> knot_vector;
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; i++) {
      knot_vector[i] = std::make_shared<baf::KnotVector>(baf::KnotVector({ParametricCoordinate(0.5)}));
    }
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; i++) {
      pugi::xml_node child = spline->child("kntVecs").first_child();
      for (int j = 0; j < i; j++) {
        child = child.next_sibling();
      }
      knot_vector[i] = std::make_shared<baf::KnotVector>(
          baf::KnotVector(util::string_operations::ConvertStringVectorToNumberVector<ParametricCoordinate>(
              util::string_operations::SplitStringAtDelimiter(child.first_child().value(), ' '))));
    }
    return knot_vector;
  }

  static std::array<Degree, PARAMETRIC_DIMENSIONALITY> StringVectorToDegreeArray(
      const std::vector<std::string> &string_vector) {
    std::array<Degree, PARAMETRIC_DIMENSIONALITY> converted{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; i++) {
      auto a = string_vector[i];
      converted[i] = Degree{std::stoi(string_vector[i])};
    }
    return converted;
  }
};
}  // namespace splinelib::src::io

#endif  // SRC_IO_XML_READER_UTILS_H_
