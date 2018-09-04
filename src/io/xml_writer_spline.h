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
  explicit XMLWriterSpline(int number_of_splines) : number_of_splines_(number_of_splines) {}

  void WriteXMLFile(const char *filename) {
    pugi::xml_document doc;
    pugi::xml_node spline_list = doc.append_child("SplineList");
    spline_list.append_attribute("NumberOfSplines") = std::to_string(number_of_splines_).c_str();
    for (int i = 0; i < number_of_splines_; i++) {
      AddSpline(&spline_list, i);
    }
    doc.save_file(filename, "  ", pugi::format_indent_attributes, pugi::encoding_utf8);
  }

 protected:
  virtual double GetDegree(int spline, int dimension) = 0;
  virtual baf::KnotVector GetKnotVector(int spline, int dimension) = 0;
  virtual char GetNumberOfControlPoints(int spline) = 0;
  virtual char GetSpaceDimension(int spline) = 0;
  virtual std::array<int, DIM> GetNumberOfPointsInEachDirection(int spline) = 0;
  virtual double GetControlPoint(int spline, std::array<int, DIM>, int dimension) = 0;

  virtual void AddWeights(pugi::xml_node *spline, int number) {}

  int number_of_splines_;

 private:
  void AddSplineAttributes(pugi::xml_node *spline, int number) {
    spline->append_attribute("splDim") = DIM;
    spline->append_attribute("spaceDim") = GetSpaceDimension(number);
    spline->append_attribute("numOfCntrlPntVars") = GetSpaceDimension(number);
    spline->append_attribute("numCntrlPnts") = GetNumberOfControlPoints(number);
  }

  void AddSpline(pugi::xml_node *spline_list, int number) {
    pugi::xml_node spline = spline_list->append_child("SplineEntry");
    AddSplineAttributes(&spline, number);
    AddControlPointVarNames(&spline, number);
    AddControlPointVars(&spline, number);
    AddWeights(&spline, number);
    AddDegrees(&spline, number);
    AddKnotVectors(&spline, number);
  }

  void AddControlPointVarNames(pugi::xml_node *spline, int number) {
    pugi::xml_node names = spline->append_child("cntrlPntVarNames");
    if (GetSpaceDimension(number) == 1) {
      names.append_child(pugi::node_pcdata).set_value("\n      x\n    ");
    } else if (GetSpaceDimension(number) == 2) {
      names.append_child(pugi::node_pcdata).set_value("\n      x y\n    ");
    } else if (GetSpaceDimension(number) == 3) {
      names.append_child(pugi::node_pcdata).set_value("\n      x y z\n    ");
    }
  }

  void AddControlPointVars(pugi::xml_node *spline, int number) {
    pugi::xml_node values = spline->append_child("cntrlPntVars");
    std::string string;
    util::MultiIndexHandler<DIM> point_handler(GetNumberOfPointsInEachDirection(number));
    for (int i = 0; i < point_handler.Get1DLength(); ++i, point_handler++) {
      auto indices = point_handler.GetIndices();
      string += "\n      ";
      for (int j = 0; j < GetSpaceDimension(number); j++) {
        string += std::to_string(GetControlPoint(number, indices, j)) + "  ";
      }
    }
    values.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  void AddDegrees(pugi::xml_node *spline, int number) {
    pugi::xml_node degrees = spline->append_child("deg");
    std::string string;
    for (int i = 0; i < DIM; i++) {
      string = string + "\n      " + std::to_string(GetDegree(number, i));
    }
    degrees.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  void AddKnotVectors(pugi::xml_node *spline, int number) {
    pugi::xml_node knot_vectors = spline->append_child("kntVecs");
    for (int i = 0; i < DIM; i++) {
      pugi::xml_node knots = knot_vectors.append_child("kntVec");
      baf::KnotVector knot_vector = GetKnotVector(number, i);
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
