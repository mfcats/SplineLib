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

#include "pugixml.hpp"

#include "b_spline.h"
#include "nurbs.h"
#include "string_operations.h"

namespace io {
template<int DIM>
class XMLReaderUtils {
 public:
  static std::array<Degree, DIM> GetDegrees(pugi::xml_node *spline) {
    return StringVectorToDegreeArray(util::StringOperations::split(spline->child("deg").first_child().value(), ' '));
  }

  static KnotVectors<DIM> GetKnotVectors(pugi::xml_node *spline) {
    KnotVectors<DIM> knot_vector;
    for (int i = 0; i < DIM; i++) {
      knot_vector[i] = std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord(0.5)}));
    }
    for (int i = 0; i < DIM; i++) {
      pugi::xml_node child = spline->child("kntVecs").first_child();
      for (int j = 0; j < i; j++) {
        child = child.next_sibling();
      }
      knot_vector[i] = std::make_shared<baf::KnotVector>(
          baf::KnotVector(util::StringOperations::StringVectorToNumberVector<ParamCoord>(
              util::StringOperations::split(child.first_child().value(), ' '))));
    }
    return knot_vector;
  }

  static std::array<Degree, DIM> StringVectorToDegreeArray(const std::vector<std::string> &string_vector) {
    std::array<Degree, DIM> converted;
    for (int i = 0; i < DIM; i++) {
      auto a = string_vector[i];
      converted[i] = Degree{std::stoi(string_vector[i])};
    }
    return converted;
  }
};
}  // namespace io

#endif  // SRC_IO_XML_READER_UTILS_H_
