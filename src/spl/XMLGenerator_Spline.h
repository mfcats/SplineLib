/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_XMLGENERATOR_SPLINE_H
#define SRC_SPL_XMLGENERATOR_SPLINE_H

#include <string>

#include "pugixml.hpp"

#include "parameter_space.h"

namespace spl {
template<int DIM>
class XMLGenerator_Spline {
 public:
  XMLGenerator_Spline() = default;

  void WriteXMLFile(const char *filename) {
    pugi::xml_document doc;
    AddSplineList(&doc);
    doc.print(std::cout, "  ");
    doc.save_file(filename, "  ", pugi::format_indent_attributes, pugi::encoding_utf8);
  }

  void ReadXMLFile(const std::string &filename) {}

 protected:
  virtual char GetNumberOfControlPoints() = 0;
  virtual char GetSpaceDimension() = 0;

  // virtual const char *GetControlPoints() = 0;

  std::shared_ptr<ParameterSpace<DIM>> parameter_space_ptr;

 private:
  void AddSplineList(pugi::xml_document *doc) {
    pugi::xml_node spline_list = doc->append_child("SplineList");
    spline_list.append_attribute("NumberOfSplines") = "1";
    AddSpline(&spline_list);
  }

  void AddSplineAttributes(pugi::xml_node *spline) {
    spline->append_attribute("splDim") = DIM;
    spline->append_attribute("spaceDim") = GetSpaceDimension();
    spline->append_attribute("numOfCntrlPntVars") = DIM;
    spline->append_attribute("numCntrlPnts") = GetNumberOfControlPoints();
  }

  void AddSpline(pugi::xml_node *spline_list) {
    pugi::xml_node spline = spline_list->append_child("SplineEntry");
    AddSplineAttributes(&spline);
    AddControlPointVarNames(&spline);
    AddControlPointVars(&spline);
    AddDegrees(&spline);
    AddKnotVectors(&spline);
  }

  void AddControlPointVarNames(pugi::xml_node *spline) {
    pugi::xml_node names = spline->append_child("cntrlPntVarNames");
    if (GetSpaceDimension() == 2) {
      names.append_child(pugi::node_pcdata).set_value("\n      x y\n    ");
    } else if (GetSpaceDimension() == 3) {
      names.append_child(pugi::node_pcdata).set_value("\n      x y z\n    ");
    }
  }

  void AddControlPointVars(pugi::xml_node *spline) {
    pugi::xml_node values = spline->append_child("cntrlPntVars");
  }

  void AddDegrees(pugi::xml_node *spline) {
    pugi::xml_node degrees = spline->append_child("deg");
    std::string string;
    for (int i = 0; i < DIM; i++) {
      string = string + "\n      " + GetString(parameter_space_ptr->GetDegree(i), 1);
    }
    degrees.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  void AddKnotVectors(pugi::xml_node *spline) {
    pugi::xml_node knot_vectors = spline->append_child("kntVecs");
    for (int i = 0; i < DIM; i++) {
      pugi::xml_node knots = knot_vectors.append_child("kntVec");
      baf::KnotVector knot_vector = (*parameter_space_ptr).GetKnotVector(i);
      std::string string;
      for (ParamCoord knot : knot_vector) {
        string += "\n        " + GetString(knot.get(), 10);
      }
      knots.append_child(pugi::node_pcdata).text() = (string + "\n      ").c_str();
    }
  }

  template<class T>
  std::string GetString(T value, int precision) const {
    std::ostringstream out;
    out << std::setprecision(precision) << value;
    std::string string = std::to_string(value);
    return string;
  }
};
}  // namespace spl

#endif  // SRC_SPL_XMLGENERATOR_SPLINE_H
