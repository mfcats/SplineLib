/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_XML_WRITER_H_
#define SRC_IO_XML_WRITER_H_

#include <any>
#include <vector>

#include "pugixml.hpp"

#include "writer.h"

namespace splinelib::src::io {
class XMLWriter : public Writer {
 public:
  XMLWriter() = default;

  void WriteFile(const std::vector<std::any> &splines, const char *filename) const override;

 private:
  void AddSpline(pugi::xml_node *spline_list, const std::any &spline) const;

  void Add1DSpline(pugi::xml_node *spline_list, const std::any &spline) const;
  void Add2DSpline(pugi::xml_node *spline_list, const std::any &spline) const;
  void Add3DSpline(pugi::xml_node *spline_list, const std::any &spline) const;
  void Add4DSpline(pugi::xml_node *spline_list, const std::any &spline) const;

  void AddSplineAttributes(pugi::xml_node *spline_node,
                           int spline_dimension,
                           int space_dimension,
                           int control_points) const;

  void AddControlPointVarNames(pugi::xml_node *spline_node, int space_dimension) const;
};
}  // namespace splinelib::src::splinelib::src::io

#endif  // SRC_IO_XML_WRITER_H_
