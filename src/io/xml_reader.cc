/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "xml_reader.h"

#include <sstream>

#include "src/util/string_operations.h"
#include "xml_reader_utils.h"

namespace splinelib::src::io {
std::vector<std::any> XMLReader::ReadFile(const char *filename) {
  std::vector<std::any> vector_of_splines;
  pugi::xml_document xml_document;
  pugi::xml_parse_result result = xml_document.load_file(filename);
  if (!result) {
    throw std::runtime_error("Input file for XML reader couldn't be parsed.");
  }
  pugi::xml_node next_spline = xml_document.child("SplineList").first_child();
  while (!next_spline.empty()) {
    AddSpline(&next_spline, &vector_of_splines);
    next_spline = next_spline.next_sibling();
  }
  return vector_of_splines;
}

void XMLReader::AddSpline(pugi::xml_node *spline, std::vector<std::any> *splines) {
  std::vector<baf::ControlPoint> control_points = GetControlPoints(spline);
  int dimension = std::stoi(spline->attribute("splDim").value());
  if (dimension == 1) {
    splines->push_back(Get1DSpline(spline, control_points));
  } else if (dimension == 2) {
    splines->push_back(Get2DSpline(spline, control_points));
  } else if (dimension == 3) {
    splines->push_back(Get3DSpline(spline, control_points));
  } else if (dimension == 4) {
    splines->push_back(Get4DSpline(spline, control_points));
  }
}

std::any XMLReader::Get1DSpline(pugi::xml_node *spline, const std::vector<baf::ControlPoint> &control_points) {
  baf::KnotVectors<1> knot_vectors = XMLReaderUtils<1>::GetKnotVectors(spline);
  std::array<Degree, 1> degrees = XMLReaderUtils<1>::GetDegrees(spline);
  if (spline->child("wght").empty()) {
    return std::make_any<std::shared_ptr<spl::BSpline<1>>>(
        std::make_shared<spl::BSpline<1>>(knot_vectors, degrees, control_points));
  }
  return std::make_any<std::shared_ptr<spl::NURBS<1>>>(
      std::make_shared<spl::NURBS<1>>(knot_vectors, degrees, control_points, GetWeights(spline)));
}

std::any XMLReader::Get2DSpline(pugi::xml_node *spline, const std::vector<baf::ControlPoint> &control_points) {
  baf::KnotVectors<2> knot_vectors = XMLReaderUtils<2>::GetKnotVectors(spline);
  std::array<Degree, 2> degrees = XMLReaderUtils<2>::GetDegrees(spline);
  if (spline->child("wght").empty()) {
    return std::make_any<std::shared_ptr<spl::BSpline<2>>>(
        std::make_shared<spl::BSpline<2>>(knot_vectors, degrees, control_points));
  }
  return std::make_any<std::shared_ptr<spl::NURBS<2>>>(
      std::make_shared<spl::NURBS<2>>(knot_vectors, degrees, control_points, GetWeights(spline)));
}

std::any XMLReader::Get3DSpline(pugi::xml_node *spline, const std::vector<baf::ControlPoint> &control_points) {
  baf::KnotVectors<3> knot_vectors = XMLReaderUtils<3>::GetKnotVectors(spline);
  std::array<Degree, 3> degrees = XMLReaderUtils<3>::GetDegrees(spline);
  if (spline->child("wght").empty()) {
    return std::make_any<std::shared_ptr<spl::BSpline<3>>>(
        std::make_shared<spl::BSpline<3>>(knot_vectors, degrees, control_points));
  }
  return std::make_any<std::shared_ptr<spl::NURBS<3>>>(
      std::make_shared<spl::NURBS<3>>(knot_vectors, degrees, control_points, GetWeights(spline)));
}

std::any XMLReader::Get4DSpline(pugi::xml_node *spline, const std::vector<baf::ControlPoint> &control_points) {
  baf::KnotVectors<4> knot_vectors = XMLReaderUtils<4>::GetKnotVectors(spline);
  std::array<Degree, 4> degrees = XMLReaderUtils<4>::GetDegrees(spline);
  if (spline->child("wght").empty()) {
    return std::make_any<std::shared_ptr<spl::BSpline<4>>>(
        std::make_shared<spl::BSpline<4>>(knot_vectors, degrees, control_points));
  }
  return std::make_any<std::shared_ptr<spl::NURBS<4>>>(
      std::make_shared<spl::NURBS<4>>(knot_vectors, degrees, control_points, GetWeights(spline)));
}

std::vector<baf::ControlPoint> XMLReader::GetControlPoints(pugi::xml_node *spline) {
  std::vector<double> vars = util::string_operations::ConvertStringVectorToNumberVector<double>(
      util::string_operations::SplitStringAtDelimiter(spline->child("cntrlPntVars").first_child().value(), ' '));
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

std::vector<double> XMLReader::GetWeights(pugi::xml_node *spline) {
  return util::string_operations::ConvertStringVectorToNumberVector<double>(
      util::string_operations::SplitStringAtDelimiter(spline->child("wght").first_child().value(), ' '));
}

int XMLReader::FindCoordinatePosition(const std::string &string) {
  std::vector<std::string> vars = util::string_operations::SplitStringAtDelimiter(string, ' ');
  for (int i = 0; i < static_cast<int>(vars.size()); i++) {
    if (vars[i] == "x") {
      return i;
    }
  }
  return 0;
}
}  // namespace splinelib::src::io
