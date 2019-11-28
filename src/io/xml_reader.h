/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_IO_XML_READER_H_
#define SRC_IO_XML_READER_H_

#include <any>
#include <string>
#include <vector>

#include "external/pugixml/pugixml.hpp"

#include "src/spl/b_spline.h"
#include "src/spl/nurbs.h"
#include "src/io/reader.h"

namespace splinelib::src::io {
class XMLReader : public Reader {
 public:
  XMLReader() = default;

  std::vector<std::any> ReadFile(const char *filename) override;

 private:
  void AddSpline(pugi::xml_node *spline, std::vector<std::any> *splines);

  std::any Get1DSpline(pugi::xml_node *spline, const std::vector<spl::ControlPoint> &control_points);
  std::any Get2DSpline(pugi::xml_node *spline, const std::vector<spl::ControlPoint> &control_points);
  std::any Get3DSpline(pugi::xml_node *spline, const std::vector<spl::ControlPoint> &control_points);
  std::any Get4DSpline(pugi::xml_node *spline, const std::vector<spl::ControlPoint> &control_points);

  std::vector<spl::ControlPoint> GetControlPoints(pugi::xml_node *spline);

  std::vector<double> GetWeights(pugi::xml_node *spline);

  int FindCoordinatePosition(const std::string &string);
};
}  // namespace splinelib::src::io

#endif  // SRC_IO_XML_READER_H_
