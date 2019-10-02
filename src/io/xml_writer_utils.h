/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_XML_WRITER_UTILS_H_
#define SRC_IO_XML_WRITER_UTILS_H_

#include <any>
#include <string>

#include "pugixml.hpp"

#include "b_spline.h"
#include "nurbs.h"

namespace splinelib::src::io {
template<int DIM>
class XMLWriterUtils {
 public:
  static void AddDegrees(pugi::xml_node *spline_node, std::shared_ptr<spl::Spline<DIM>> spline_ptr) {
    pugi::xml_node degrees = spline_node->append_child("deg");
    std::string string;
    for (int i = 0; i < DIM; i++) {
      string = string + "\n      " + std::to_string(spline_ptr->GetDegree(i).get());
    }
    degrees.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  static void AddKnotVectors(pugi::xml_node *spline_node, std::shared_ptr<spl::Spline<DIM>> spline_ptr) {
    pugi::xml_node knot_vectors = spline_node->append_child("kntVecs");
    for (int i = 0; i < DIM; i++) {
      pugi::xml_node knots = knot_vectors.append_child("kntVec");
      baf::KnotVector knot_vector = *spline_ptr->GetKnotVector(i);
      std::string string;
      for (ParametricCoordinate knot : knot_vector) {
        string += "\n        " + std::to_string(knot.get());
      }
      knots.append_child(pugi::node_pcdata).text() = (string + "\n      ").c_str();
    }
  }

  static void AddControlPointVars(pugi::xml_node *spline_node, std::shared_ptr<spl::Spline<DIM>> spline_ptr) {
    pugi::xml_node values = spline_node->append_child("cntrlPntVars");
    std::string string;
    util::MultiIndexHandler<DIM> point_handler(spline_ptr->GetPointsPerDirection());
    for (int i = 0; i < point_handler.Get1DLength(); ++i, point_handler++) {
      auto indices = point_handler.GetIndices();
      string += "\n      ";
      for (int j = 0; j < spline_ptr->GetPointDim(); j++) {
        string += std::to_string(spline_ptr->GetControlPoint(indices, j)) + "  ";
      }
    }
    values.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  static void AddWeights(pugi::xml_node *spline_node, const std::any &spline) {
    pugi::xml_node weights = spline_node->append_child("wght");
    std::string string;
    std::shared_ptr<spl::NURBS<DIM>> nurbs = std::any_cast<std::shared_ptr<spl::NURBS<DIM>>>(spline);
    util::MultiIndexHandler<DIM> weight_handler(nurbs->GetPointsPerDirection());
    for (int i = 0; i < weight_handler.Get1DLength(); ++i, weight_handler++) {
      auto indices = weight_handler.GetIndices();
      string += "\n      " + std::to_string(nurbs->GetWeight(indices)) + "  ";
    }
    weights.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }
};
}  // namespace splinelib::src::io

#endif  // SRC_IO_XML_WRITER_UTILS_H_
