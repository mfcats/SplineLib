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
    doc.document_element().append_attribute("encoding") = "utf-8";
    pugi::xml_node list = doc.append_child("SplineList");
    list.append_attribute("NumberOfSplines") = "1";
    pugi::xml_node spline = list.append_child("SplineEntry");
    spline.append_attribute("splDim") = DIM;
    spline.append_attribute("spaceDim") = GetSpaceDimension();
    spline.append_attribute("numOfCntrlPntVars") = DIM;
    spline.append_attribute("numCntrlPnts") = GetNumberOfControlPoints();
    pugi::xml_node names = spline.append_child("cntrlPntVarNames");
    names.append_child(pugi::node_pcdata).set_value("\n      x y z\n    ");
    pugi::xml_node values = spline.append_child("cntrlPntVars");
    pugi::xml_node degrees = spline.append_child("deg");
    degrees.append_child(pugi::node_pcdata).text() = GetDegrees();
    pugi::xml_node knot_vectors = spline.append_child("kntVecs");
    for (int i = 0; i < DIM; i++) {
      pugi::xml_node knots = knot_vectors.append_child("kntVec");
      baf::KnotVector knot_vector = (*parameter_space_ptr).GetKnotVector(i);
      knots.append_child(pugi::node_pcdata).text() = ConvertKnotVectorToString(knot_vector);
    }
    doc.print(std::cout, "  ");
    std::cout << "Saving result: " << doc.save_file(filename, "  ", pugi::format_indent_attributes, pugi::encoding_utf8)
              << std::endl;
  }

  void ReadXMLFile(const std::string &filename) {}

 protected:
  virtual char GetNumberOfControlPoints() = 0;
  virtual char GetSpaceDimension() = 0;

  const char *GetDegrees() {
    std::string string = "test";
    for (int i = 0; i < DIM; i++) {
      //string = string + "\n      " + GetString(parameter_space_ptr->GetDegree(i), 1);
    }
    // string += "test"; //"\n      ";
    return string.c_str();
  }

  virtual const char *GetControlPoints() = 0;
  const char *ConvertKnotVectorToString(baf::KnotVector knot_vector) {
    std::string string;
    for (ParamCoord knot : knot_vector) {
      string = string + "\n        " + GetString(knot.get(), 10);
    }
    string += "\n       ";
    return string.c_str();
  }

  template<class T>
  std::string GetString(T value, int precision) const {
    std::ostringstream out;
    out << std::setprecision(precision) << value;
    std::string string = std::to_string(value);
    return string;
  }

  std::shared_ptr<ParameterSpace<DIM>> parameter_space_ptr;
};
}  // namespace spl

#endif  // SRC_SPL_XMLGENERATOR_SPLINE_H
