/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_XML_READER_H_
#define SRC_IO_XML_READER_H_

#include <any>
#include <sstream>
#include <string>
#include <vector>

#include "pugixml.hpp"

#include "b_spline.h"
#include "nurbs.h"

namespace io {
template<int DIM>
class XMLReader {
 public:
  XMLReader() = default;

  std::any ReadXMLFile(const std::string &filename) {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());
    if (!result) {
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
    spl::ParameterSpace<DIM> parameterSpace = spl::ParameterSpace<DIM>(knot_vector, degree);
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
    if (spline.child("wght").empty()) {
      return std::make_any<spl::BSpline<DIM>>(parameterSpace, physical_space);
    } else {
      std::vector<double> weights = StringVectorToDoubleVector(split(spline.child("wght").first_child().value()));
      spl::WeightedPhysicalSpace<DIM> weightedPhysicalSpace(control_points_, weights, number_of_control_points);
      return std::make_any<spl::NURBS<DIM>>(parameterSpace, weightedPhysicalSpace);
    }
  }

 private:
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
      converted.push_back(strtod(string.c_str(), nullptr));
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
      converted.emplace_back(strtod(string.c_str(), nullptr));
    }
    return baf::KnotVector(converted);
  }

  int FindCoordinatePosition(const std::string &string) {
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
}  // namespace io

#endif  // SRC_IO_XML_READER_H_
