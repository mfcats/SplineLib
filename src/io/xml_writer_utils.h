/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_IO_XML_WRITER_UTILS_H_
#define SRC_IO_XML_WRITER_UTILS_H_

#include <any>
#include <string>

#include "external/pugixml/pugixml.hpp"

#include "src/spl/b_spline.h"
#include "src/spl/nurbs.h"
#include "src/util/string_operations.h"

namespace splinelib::src::io {
template<int PARAMETRIC_DIMENSIONALITY>
class XMLWriterUtils {
 public:
  static void AddDegrees(pugi::xml_node *spline_node,
      std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>> spline_ptr) {
    pugi::xml_node degrees = spline_node->append_child("deg");
    std::string string;
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; i++) {
      string += "\n      " + std::to_string(spline_ptr->GetDegree(i).Get());
    }
    degrees.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  static void AddKnotVectors(pugi::xml_node *spline_node,
      std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>> spline_ptr) {
    pugi::xml_node knot_vectors = spline_node->append_child("kntVecs");
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; i++) {
      pugi::xml_node knots = knot_vectors.append_child("kntVec");
      baf::KnotVector knot_vector = *spline_ptr->GetKnotVector(i);
      std::string string;
      for (ParametricCoordinate knot : knot_vector) {
        string += "\n        " + util::string_operations::GetStringWithHighPrecision(knot.Get());
      }
      knots.append_child(pugi::node_pcdata).text() = (string + "\n      ").c_str();
    }
  }

  static void AddControlPointVars(pugi::xml_node *spline_node,
      std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>> spline_ptr) {
    pugi::xml_node values = spline_node->append_child("cntrlPntVars");
    std::string string;
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(spline_ptr->GetNumberOfPointsPerDirection());
    for (int i = 0; i < point_handler.GetNumberOfTotalMultiIndices(); ++i, point_handler++) {
      auto indices = point_handler.GetCurrentIndex();
      string += "\n      ";
      for (int j = 0; j < spline_ptr->GetPointDim(); j++) {
        string += util::string_operations::GetStringWithHighPrecision(spline_ptr->GetControlPoint(indices, j)) + "  ";
      }
    }
    values.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  static void AddWeights(pugi::xml_node *spline_node, const std::any &spline) {
    pugi::xml_node weights = spline_node->append_child("wght");
    std::string string;
    auto nurbs = std::any_cast<std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>>>(spline);
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> weight_handler(nurbs->GetNumberOfPointsPerDirection());
    for (int i = 0; i < weight_handler.GetNumberOfTotalMultiIndices(); ++i, weight_handler++) {
      auto indices = weight_handler.GetCurrentIndex();
      string += "\n      " + util::string_operations::GetStringWithHighPrecision(nurbs->GetWeight(indices)) + "  ";
    }
    weights.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }
};
}  // namespace splinelib::src::io

#endif  // SRC_IO_XML_WRITER_UTILS_H_
