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
#include "string_operations.h"

namespace io {
template<int DIM>
class XMLReader {
 public:
  XMLReader() = default;

  std::vector<std::any> ReadXMLFile(const std::string &filename) {
    std::vector<std::any> vector_of_splines;
    pugi::xml_document xml_document;
    pugi::xml_parse_result result = xml_document.load_file(filename.c_str());
    if (!result) {
      throw std::runtime_error("File couldn't be loaded.");
    }
    pugi::xml_node next_spline = xml_document.child("SplineList").first_child();
    while (!next_spline.empty()) {
      AddSpline(&next_spline, &vector_of_splines);
      next_spline = next_spline.next_sibling();
    }
    return vector_of_splines;
  }

 private:
  void AddSpline(pugi::xml_node *spline, std::vector<std::any> *splines) {
    spl::ParameterSpace<DIM> parameter_space = GetParameterSpace(spline);
    std::vector<baf::ControlPoint> control_points_ = GetControlPoints(spline);
    std::array<int, DIM> number_of_control_points = GetNumberOfControlPoints(parameter_space);
    if (spline->child("wght").empty()) {
      splines->push_back(std::make_any<spl::BSpline<DIM>>(parameter_space, spl::PhysicalSpace<DIM>(
          control_points_, number_of_control_points)));
    } else {
      std::vector<double> weights = util::StringOperations::StringVectorToNumberVector<double>(
          util::StringOperations::split(spline->child("wght").first_child().value(), ' '));
      spl::WeightedPhysicalSpace<DIM> weightedPhysicalSpace(control_points_, weights, number_of_control_points);
      splines->push_back(std::make_any<spl::NURBS<DIM>>(parameter_space, weightedPhysicalSpace));
    }
  }

  std::vector<baf::ControlPoint> GetControlPoints(pugi::xml_node *spline) {
    std::vector<double> vars =
        util::StringOperations::StringVectorToNumberVector<double>(util::StringOperations::split(spline->child(
            "cntrlPntVars").first_child().value(), ' '));
    int start = FindCoordinatePosition(spline->child("cntrlPntVarNames").first_child().value());
    int dimension = std::stoi(spline->attribute("spaceDim").value());
    int numberOfVars = std::stoi(spline->attribute("numOfCntrlPntVars").value());
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

  std::array<int, DIM> GetNumberOfControlPoints(spl::ParameterSpace<DIM> parameter_space) {
    std::array<int, DIM> number_of_control_points;
    for (int i = 0; i < DIM; i++) {
      auto h = parameter_space.GetKnotVector(i).GetNumberOfKnots() - parameter_space.GetDegree(i) - 1;
      number_of_control_points[i] =
          parameter_space.GetKnotVector(i).GetNumberOfKnots() - parameter_space.GetDegree(i) - 1;
    }
    return number_of_control_points;
  };

  spl::ParameterSpace<DIM> GetParameterSpace(pugi::xml_node *spline) {
    std::array<int, DIM>
        degree = StringVectorToIntArray(util::StringOperations::split(spline->child("deg").first_child().value(), ' '));
    std::array<baf::KnotVector, DIM> knot_vector;
    for (int i = 0; i < DIM; i++) {
      knot_vector[i] = baf::KnotVector({ParamCoord(0.5)});
    }
    for (int i = 0; i < DIM; i++) {
      pugi::xml_node child = spline->child("kntVecs").first_child();
      for (int j = 0; j < i; j++) {
        child = child.next_sibling();
      }
      knot_vector[i] = baf::KnotVector(util::StringOperations::StringVectorToNumberVector<ParamCoord>(
          util::StringOperations::split(child.first_child().value(), ' ')));
    }
    return spl::ParameterSpace<DIM>(knot_vector, degree);
  }

  std::array<int, DIM> StringVectorToIntArray(const std::vector<std::string> &string_vector) {
    std::array<int, DIM> converted;
    for (int i = 0; i < DIM; i++) {
      converted[i] = std::stoi(string_vector[i]);
    }
    return converted;
  }

  int FindCoordinatePosition(const std::string &string) {
    std::vector<std::string> vars = util::StringOperations::split(string, ' ');
    for (int i = 0; i < vars.size(); i++) {
      if (vars[i] == "x") {
        return i;
      }
    }
    return 0;
  }
};
}  // namespace io

#endif  // SRC_IO_XML_READER_H_
