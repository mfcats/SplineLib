/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_XML_WRITER_SPLINE_H_
#define SRC_IO_XML_WRITER_SPLINE_H_

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

#include "pugixml.hpp"

#include "b_spline.h"
#include "nurbs.h"

namespace io {
template<int DIM>
class XMLWriterSpline {
 public:
  XMLWriterSpline() = default;

  void WriteXMLFile(const char *filename) {
    pugi::xml_document doc;
    AddSplineList(&doc);
    doc.save_file(filename, "  ", pugi::format_indent_attributes, pugi::encoding_utf8);
  }

 protected:
  virtual double GetDegree(int dimension) = 0;
  virtual baf::KnotVector GetKnotVector(int dimension) = 0;
  virtual char GetNumberOfControlPoints() = 0;
  virtual char GetSpaceDimension() = 0;
  virtual std::array<int, DIM> GetNumberOfPointsInEachDirection() = 0;
  virtual double GetControlPoint(std::array<int, DIM>, int dimension) = 0;

  virtual void AddWeights(pugi::xml_node *spline) {}

 private:
  void AddSplineList(pugi::xml_document *doc) {
    pugi::xml_node spline_list = doc->append_child("SplineList");
    spline_list.append_attribute("NumberOfSplines") = "1";
    AddSpline(&spline_list);
  }

  void AddSplineAttributes(pugi::xml_node *spline) {
    spline->append_attribute("splDim") = DIM;
    spline->append_attribute("spaceDim") = GetSpaceDimension();
    spline->append_attribute("numOfCntrlPntVars") = GetSpaceDimension();
    spline->append_attribute("numCntrlPnts") = GetNumberOfControlPoints();
  }

  void AddSpline(pugi::xml_node *spline_list) {
    pugi::xml_node spline = spline_list->append_child("SplineEntry");
    AddSplineAttributes(&spline);
    AddControlPointVarNames(&spline);
    AddControlPointVars(&spline);
    AddWeights(&spline);
    AddDegrees(&spline);
    AddKnotVectors(&spline);
  }

  void AddControlPointVarNames(pugi::xml_node *spline) {
    pugi::xml_node names = spline->append_child("cntrlPntVarNames");
    if (GetSpaceDimension() == 1) {
      names.append_child(pugi::node_pcdata).set_value("\n      x\n    ");
    } else if (GetSpaceDimension() == 2) {
      names.append_child(pugi::node_pcdata).set_value("\n      x y\n    ");
    } else if (GetSpaceDimension() == 3) {
      names.append_child(pugi::node_pcdata).set_value("\n      x y z\n    ");
    }
  }

  void AddControlPointVars(pugi::xml_node *spline) {
    pugi::xml_node values = spline->append_child("cntrlPntVars");
    std::string string;
    util::MultiIndexHandler<DIM> point_handler(GetNumberOfPointsInEachDirection());
    for (int i = 0; i < point_handler.Get1DLength(); ++i, point_handler++) {
      auto indices = point_handler.GetIndices();
      string += "\n      ";
      for (int j = 0; j < GetSpaceDimension(); j++) {
        string += std::to_string(GetControlPoint(indices, j)) + "  ";
      }
    }
    values.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  void AddDegrees(pugi::xml_node *spline) {
    pugi::xml_node degrees = spline->append_child("deg");
    std::string string;
    for (int i = 0; i < DIM; i++) {
      string = string + "\n      " + std::to_string(GetDegree(i));
    }
    degrees.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  void AddKnotVectors(pugi::xml_node *spline) {
    pugi::xml_node knot_vectors = spline->append_child("kntVecs");
    for (int i = 0; i < DIM; i++) {
      pugi::xml_node knots = knot_vectors.append_child("kntVec");
      baf::KnotVector knot_vector = GetKnotVector(i);
      std::string string;
      for (ParamCoord knot : knot_vector) {
        string += "\n        " + std::to_string(knot.get());
      }
      knots.append_child(pugi::node_pcdata).text() = (string + "\n      ").c_str();
    }
  }
};
}  // namespace io

#endif  // SRC_IO_XML_WRITER_SPLINE_H_
