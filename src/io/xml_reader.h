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

  std::vector<std::any> ReadXMLFile(const char *filename) {
    std::vector<std::any> vector_of_splines;
    pugi::xml_document xml_document;
    pugi::xml_parse_result result = xml_document.load_file(filename);
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
    std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vectors = GetKnotVectors(spline);
    std::array<Degree, DIM> degrees = GetDegrees(spline);
    std::vector<baf::ControlPoint> control_points = GetControlPoints(spline);
    if (spline->child("wght").empty()) {
      splines->push_back(std::make_any<std::shared_ptr<spl::BSpline<DIM>>>(
          std::make_shared<spl::BSpline<DIM>>(knot_vectors, degrees, control_points)));
    } else {
      splines->push_back(std::make_any<std::shared_ptr<spl::NURBS<DIM>>>(
          std::make_shared<spl::NURBS<DIM>>(knot_vectors, degrees, control_points, GetWeights(spline))));
    }
  }

  std::vector<baf::ControlPoint> GetControlPoints(pugi::xml_node *spline) {
    std::vector<double> vars = util::StringOperations::StringVectorToNumberVector<double>(
        util::StringOperations::split(spline->child("cntrlPntVars").first_child().value(), ' '));
    int start = FindCoordinatePosition(spline->child("cntrlPntVarNames").first_child().value());
    int dimension = std::stoi(spline->attribute("spaceDim").value());
    int number_of_vars = std::stoi(spline->attribute("numOfCntrlPntVars").value());
    int number_of_points = std::stoi(spline->attribute("numCntrlPnts").value());
    std::vector<baf::ControlPoint> points;
    for (int i = 0; i < number_of_points; i++) {
      std::vector<double> coordinates;
      for (int j = start; j < start + dimension; j++) {
        coordinates.push_back(vars[i * number_of_vars + j]);
      }
      points.emplace_back(coordinates);
    }
    return points;
  }

  std::array<Degree, DIM> GetDegrees(pugi::xml_node *spline) {
    return StringVectorToDegreeArray(util::StringOperations::split(spline->child("deg").first_child().value(), ' '));
  }

  std::array<std::shared_ptr<baf::KnotVector>, DIM> GetKnotVectors(pugi::xml_node *spline) {
    std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vector;
    for (int i = 0; i < DIM; i++) {
      knot_vector[i] = std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord(0.5)}));
    }
    for (int i = 0; i < DIM; i++) {
      pugi::xml_node child = spline->child("kntVecs").first_child();
      for (int j = 0; j < i; j++) {
        child = child.next_sibling();
      }
      knot_vector[i] = std::make_shared<baf::KnotVector>(
          baf::KnotVector(util::StringOperations::StringVectorToNumberVector<ParamCoord>(
              util::StringOperations::split(child.first_child().value(), ' '))));
    }
    return knot_vector;
  }

  std::vector<double> GetWeights(pugi::xml_node *spline) {
    return util::StringOperations::StringVectorToNumberVector<double>(
        util::StringOperations::split(spline->child("wght").first_child().value(), ' '));
  }

  std::array<Degree, DIM> StringVectorToDegreeArray(const std::vector<std::string> &string_vector) {
    std::array<Degree, DIM> converted;
    for (int i = 0; i < DIM; i++) {
      converted[i] = Degree{std::stoi(string_vector[i])};
    }
    return converted;
  }

  int FindCoordinatePosition(const std::string &string) {
    std::vector<std::string> vars = util::StringOperations::split(string, ' ');
    for (int i = 0; i < static_cast<int>(vars.size()); i++) {
      if (vars[i] == "x") {
        return i;
      }
    }
    return 0;
  }
};
}  // namespace io

#endif  // SRC_IO_XML_READER_H_
