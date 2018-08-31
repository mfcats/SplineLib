/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_XML_GENERATOR_SPLINE_H_
#define SRC_SPL_XML_GENERATOR_SPLINE_H_

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

#include "pugixml.hpp"

#include "b_spline.h"

namespace spl {
template<int DIM>
class XMLGenerator_Spline {
 public:
  XMLGenerator_Spline() = default;

  void WriteXMLFile(const char *filename) {
    pugi::xml_document doc;
    AddSplineList(&doc);
    doc.save_file(filename, "  ", pugi::format_indent_attributes, pugi::encoding_utf8);
  }

  std::shared_ptr<BSpline<DIM>> ReadXMLFile(const char *filename) {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename);
    if (!result) {
      std::cout << result.description();
      throw std::runtime_error("File couldn't be loaded.");
    }
    pugi::xml_node spline = doc.child("SplineList").first_child();
    std::array<int, DIM> degree = StringVectorToIntArray(split(spline.child("deg").first_child().value()));
    std::array<baf::KnotVector, DIM> knot_vector;
    for (int i = 0; i < DIM; i++) {
      knot_vector[i] = baf::KnotVector({ParamCoord(0.5)});
    }
    for (int i = 0; i < DIM; i++) {
      pugi::xml_node child = spline.child("kntVecs").first_child();
      for (int j = 0; j < i; j++) {
        child = child.next_sibling();
      }
      knot_vector[i] = StringVectorToKnotVector(split(child.first_child().value()));
    }
    ParameterSpace<DIM> parameterSpace = ParameterSpace<DIM>(knot_vector, degree);
    // std::vector<double> weights;
    int dim = std::stoi(spline.attribute("spaceDim").value());
    std::vector<double>
        control_point_vars = StringVectorToDoubleVector(split(spline.child("cntrlPntVars").first_child().value()));
    std::string var_names = spline.child("cntrlPntVarNames").first_child().value();
    int dimension = std::stoi(spline.attribute("spaceDim").value());
    int numberOfVars = std::stoi(spline.attribute("numOfCntrlPntVars").value());
    std::vector<baf::ControlPoint> control_points_ =
        GetControlPoints(control_point_vars, FindCoordinatePosition(var_names), dimension, numberOfVars);
    std::array<int, DIM> number_of_control_points;
    for (int i = 0; i < DIM; i++) {
      auto h = knot_vector[i].GetNumberOfKnots() - degree[i] - 1;
      number_of_control_points[i] = knot_vector[i].GetNumberOfKnots() - degree[i] - 1;
    }
    spl::PhysicalSpace<DIM> physical_space = spl::PhysicalSpace<DIM>(control_points_, number_of_control_points);
    std::shared_ptr<BSpline<DIM>> b_spline = std::make_shared<BSpline<DIM>>(parameterSpace, physical_space);
    return b_spline;
  }

 protected:
  virtual char GetNumberOfControlPoints() = 0;
  virtual char GetSpaceDimension() = 0;
  virtual std::array<int, DIM> GetNumberOfPointsInEachDirection() = 0;
  virtual double GetControlPoint(std::array<int, DIM>, int dimension) = 0;

  virtual void AddWeights(pugi::xml_node *spline) {}

  template<class T>
  std::string GetString(T value) const {
    std::string string = std::to_string(value);
    return string;
  }

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
    if (GetSpaceDimension() == 2) {
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
        string += GetString(GetControlPoint(indices, j)) + "  ";
      }
    }
    values.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  void AddDegrees(pugi::xml_node *spline) {
    pugi::xml_node degrees = spline->append_child("deg");
    std::string string;
    for (int i = 0; i < DIM; i++) {
      string = string + "\n      " + GetString(parameter_space_ptr->GetDegree(i));
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
        string += "\n        " + GetString(knot.get());
      }
      knots.append_child(pugi::node_pcdata).text() = (string + "\n      ").c_str();
    }
  }

  std::vector<std::string> split(const std::string &string) {
    std::stringstream ss(string);
    std::string line;
    std::string value;
    std::vector<std::string> splitted;
    while (std::getline(ss, line)) {
      std::stringstream test(line);
      while (std::getline(test, value, ' ')) {
        if (!value.empty()) splitted.push_back(value);
      }
    }
    return splitted;
  }

  std::vector<double> StringVectorToDoubleVector(const std::vector<std::string> &string_vector) {
    std::vector<double> converted;
    for (const std::string &string : string_vector) {
      converted.push_back(atof(string.c_str()));
    }
    return converted;
  }

  std::array<int, DIM> StringVectorToIntArray(const std::vector<std::string> &string_vector) {
    std::array<int, DIM> converted;
    for (int i = 0; i < DIM; i++) {
      converted[i] = std::stoi(string_vector[i]);
    }
    return converted;
  }

  baf::KnotVector StringVectorToKnotVector(const std::vector<std::string> &string_vector) {
    std::vector<ParamCoord> converted;
    for (const std::string &string : string_vector) {
      converted.emplace_back(atof(string.c_str()));
    }
    return baf::KnotVector(converted);
  }

  int FindCoordinatePosition(std::string string) {
    std::vector<std::string> vars = split(string);
    for (int i = 0; i < vars.size(); i++) {
      if (vars[i] == "x") {
        return i;
      }
    }
    return 0;
  }

  std::vector<baf::ControlPoint> GetControlPoints(std::vector<double> vars,
                                                  int start,
                                                  int dimension,
                                                  int numberOfVars) {
    std::vector<baf::ControlPoint> points;
    for (int i = 0; i < vars.size() / numberOfVars; i++) {
      std::vector<double> coordinates;
      for (int j = start; j < start + dimension; j++) {
        coordinates.push_back(vars[i * numberOfVars + j]);
      }
      points.emplace_back(coordinates);
    }
    return points;
  }
};
}  // namespace spl

#endif  // SRC_SPL_XML_GENERATOR_SPLINE_H_
