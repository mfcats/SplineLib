/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_XML_WRITER_H_
#define SRC_IO_XML_WRITER_H_

#include <any>
#include <string>
#include <vector>

#include "pugixml.hpp"

#include "any_casts.h"
#include "b_spline.h"
#include "nurbs.h"

namespace io {
template<int DIM>
class XMLWriter {
 public:
  explicit XMLWriter(const std::vector<std::any> &splines) : splines_(splines) {}

  void WriteXMLFile(const char *filename) {
    pugi::xml_document doc;
    pugi::xml_node spline_list = doc.append_child("SplineList");
    spline_list.append_attribute("NumberOfSplines") = std::to_string(splines_.size()).c_str();
    for (int i = 0; i < static_cast<int>(splines_.size()); i++) {
      AddSpline(&spline_list, i);
    }
    doc.save_file(filename, "  ", pugi::format_indent_attributes, pugi::encoding_utf8);
  }

 private:
  void AddSpline(pugi::xml_node *spline_list, int spline_number) {
    std::shared_ptr<spl::Spline<DIM>> spline_ptr = util::AnyCasts<DIM>::GetSpline(splines_[spline_number]);
    pugi::xml_node spline_node = spline_list->append_child("SplineEntry");
    AddSplineAttributes(&spline_node, spline_ptr);
    AddControlPointVarNames(&spline_node, spline_ptr);
    AddControlPointVars(&spline_node, spline_ptr);
    if (util::AnyCasts<DIM>::IsRational(splines_[spline_number])) AddWeights(&spline_node, spline_number);
    AddDegrees(&spline_node, spline_ptr);
    AddKnotVectors(&spline_node, spline_ptr);
  }

  void AddSplineAttributes(pugi::xml_node *spline_node, std::shared_ptr<spl::Spline<DIM>> spline_ptr) {
    spline_node->append_attribute("splDim") = DIM;
    spline_node->append_attribute("spaceDim") = spline_ptr->GetDimension();
    spline_node->append_attribute("numOfCntrlPntVars") = spline_ptr->GetDimension();
    spline_node->append_attribute("numCntrlPnts") = spline_ptr->GetNumberOfControlPoints();
  }

  void AddControlPointVarNames(pugi::xml_node *spline_node, std::shared_ptr<spl::Spline<DIM>> spline_ptr) {
    pugi::xml_node names = spline_node->append_child("cntrlPntVarNames");
    if (spline_ptr->GetDimension() == 1) {
      names.append_child(pugi::node_pcdata).set_value("\n      x\n    ");
    } else if (spline_ptr->GetDimension() == 2) {
      names.append_child(pugi::node_pcdata).set_value("\n      x y\n    ");
    } else if (spline_ptr->GetDimension() == 3) {
      names.append_child(pugi::node_pcdata).set_value("\n      x y z\n    ");
    }
  }

  void AddControlPointVars(pugi::xml_node *spline_node, std::shared_ptr<spl::Spline<DIM>> spline_ptr) {
    pugi::xml_node values = spline_node->append_child("cntrlPntVars");
    std::string string;
    util::MultiIndexHandler<DIM> point_handler(spline_ptr->GetPointsPerDirection());
    for (int i = 0; i < point_handler.Get1DLength(); ++i, point_handler++) {
      auto indices = point_handler.GetIndices();
      string += "\n      ";
      for (int j = 0; j < spline_ptr->GetDimension(); j++) {
        string += std::to_string(spline_ptr->GetControlPoint(indices, j)) + "  ";
      }
    }
    values.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  void AddDegrees(pugi::xml_node *spline_node, std::shared_ptr<spl::Spline<DIM>> spline_ptr) {
    pugi::xml_node degrees = spline_node->append_child("deg");
    std::string string;
    for (int i = 0; i < DIM; i++) {
      string = string + "\n      " + std::to_string(spline_ptr->GetDegree(i).get());
    }
    degrees.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  void AddKnotVectors(pugi::xml_node *spline_node, std::shared_ptr<spl::Spline<DIM>> spline_ptr) {
    pugi::xml_node knot_vectors = spline_node->append_child("kntVecs");
    for (int i = 0; i < DIM; i++) {
      pugi::xml_node knots = knot_vectors.append_child("kntVec");
      baf::KnotVector knot_vector = *spline_ptr->GetKnotVector(i);
      std::string string;
      for (ParamCoord knot : knot_vector) {
        string += "\n        " + std::to_string(knot.get());
      }
      knots.append_child(pugi::node_pcdata).text() = (string + "\n      ").c_str();
    }
  }

  void AddWeights(pugi::xml_node *spline_node, int spline_number) {
    pugi::xml_node weights = spline_node->append_child("wght");
    std::string string;
    std::shared_ptr<spl::NURBS<DIM>> nurbs = std::any_cast<std::shared_ptr<spl::NURBS<DIM>>>(splines_[spline_number]);
    util::MultiIndexHandler<DIM> weight_handler(nurbs->GetPointsPerDirection());
    for (int i = 0; i < weight_handler.Get1DLength(); ++i, weight_handler++) {
      auto indices = weight_handler.GetIndices();
      string += "\n      " + std::to_string(nurbs->GetWeight(indices)) + "  ";
    }
    weights.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  std::vector<std::any> splines_;
};
}  // namespace io

#endif  // SRC_IO_XML_WRITER_H_
