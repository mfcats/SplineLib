/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "xml_writer.h"

#include <string>

#include "src/util/any_casts.h"
#include "b_spline.h"
#include "nurbs.h"
#include "xml_writer_utils.h"

namespace splinelib::src::io {
void XMLWriter::WriteFile(const std::vector<std::any> &splines, const char *filename) const {
  pugi::xml_document doc;
  pugi::xml_node spline_list = doc.append_child("SplineList");
  spline_list.append_attribute("NumberOfSplines") = std::to_string(splines.size()).c_str();
  for (const auto &spline : splines) {
    AddSpline(&spline_list, spline);
  }
  doc.save_file(filename, "  ", pugi::format_indent_attributes, pugi::encoding_utf8);
}

void XMLWriter::AddSpline(pugi::xml_node *spline_list, const std::any &spline) const {
  int spline_dimension = util::AnyCasts::GetSplineDimension(spline);
  if (spline_dimension == 1) {
    Add1DSpline(spline_list, spline);
  } else if (spline_dimension == 2) {
    Add2DSpline(spline_list, spline);
  } else if (spline_dimension == 3) {
    Add3DSpline(spline_list, spline);
  } else if (spline_dimension == 4) {
    Add4DSpline(spline_list, spline);
  }
}

void XMLWriter::Add1DSpline(pugi::xml_node *spline_list, const std::any &spline) const {
  std::shared_ptr<spl::Spline<1>> spline_ptr = util::AnyCasts::GetSpline<1>(spline);
  pugi::xml_node spline_node = spline_list->append_child("SplineEntry");
  AddSplineAttributes(&spline_node, 1, spline_ptr->GetPointDim(), spline_ptr->GetNumberOfControlPoints());
  AddControlPointVarNames(&spline_node, spline_ptr->GetPointDim());
  XMLWriterUtils<1>::AddControlPointVars(&spline_node, spline_ptr);
  if (util::AnyCasts::IsRational<1>(spline)) {
    XMLWriterUtils<1>::AddWeights(&spline_node, spline);
  }
  XMLWriterUtils<1>::AddDegrees(&spline_node, spline_ptr);
  XMLWriterUtils<1>::AddKnotVectors(&spline_node, spline_ptr);
}

void XMLWriter::Add2DSpline(pugi::xml_node *spline_list, const std::any &spline) const {
  std::shared_ptr<spl::Spline<2>> spline_ptr = util::AnyCasts::GetSpline<2>(spline);
  pugi::xml_node spline_node = spline_list->append_child("SplineEntry");
  AddSplineAttributes(&spline_node, 2, spline_ptr->GetPointDim(), spline_ptr->GetNumberOfControlPoints());
  AddControlPointVarNames(&spline_node, spline_ptr->GetPointDim());
  XMLWriterUtils<2>::AddControlPointVars(&spline_node, spline_ptr);
  if (util::AnyCasts::IsRational<2>(spline)) {
    XMLWriterUtils<2>::AddWeights(&spline_node, spline);
  }
  XMLWriterUtils<2>::AddDegrees(&spline_node, spline_ptr);
  XMLWriterUtils<2>::AddKnotVectors(&spline_node, spline_ptr);
}

void XMLWriter::Add3DSpline(pugi::xml_node *spline_list, const std::any &spline) const {
  std::shared_ptr<spl::Spline<3>> spline_ptr = util::AnyCasts::GetSpline<3>(spline);
  pugi::xml_node spline_node = spline_list->append_child("SplineEntry");
  AddSplineAttributes(&spline_node, 3, spline_ptr->GetPointDim(), spline_ptr->GetNumberOfControlPoints());
  AddControlPointVarNames(&spline_node, spline_ptr->GetPointDim());
  XMLWriterUtils<3>::AddControlPointVars(&spline_node, spline_ptr);
  if (util::AnyCasts::IsRational<3>(spline)) {
    XMLWriterUtils<3>::AddWeights(&spline_node, spline);
  }
  XMLWriterUtils<3>::AddDegrees(&spline_node, spline_ptr);
  XMLWriterUtils<3>::AddKnotVectors(&spline_node, spline_ptr);
}

void XMLWriter::Add4DSpline(pugi::xml_node *spline_list, const std::any &spline) const {
  std::shared_ptr<spl::Spline<4>> spline_ptr = util::AnyCasts::GetSpline<4>(spline);
  pugi::xml_node spline_node = spline_list->append_child("SplineEntry");
  AddSplineAttributes(&spline_node, 4, spline_ptr->GetPointDim(), spline_ptr->GetNumberOfControlPoints());
  AddControlPointVarNames(&spline_node, spline_ptr->GetPointDim());
  XMLWriterUtils<4>::AddControlPointVars(&spline_node, spline_ptr);
  if (util::AnyCasts::IsRational<4>(spline)) {
    XMLWriterUtils<4>::AddWeights(&spline_node, spline);
  }
  XMLWriterUtils<4>::AddDegrees(&spline_node, spline_ptr);
  XMLWriterUtils<4>::AddKnotVectors(&spline_node, spline_ptr);
}

void XMLWriter::AddSplineAttributes(pugi::xml_node *spline_node,
                                        int spline_dimension,
                                        int space_dimension,
                                        int control_points) const {
  spline_node->append_attribute("splDim") = spline_dimension;
  spline_node->append_attribute("spaceDim") = space_dimension;
  spline_node->append_attribute("numOfCntrlPntVars") = space_dimension;
  spline_node->append_attribute("numCntrlPnts") = control_points;
}

void XMLWriter::AddControlPointVarNames(pugi::xml_node *spline_node, int space_dimension) const {
  pugi::xml_node names = spline_node->append_child("cntrlPntVarNames");
  if (space_dimension == 1) {
    names.append_child(pugi::node_pcdata).set_value("\n      x\n    ");
  } else if (space_dimension == 2) {
    names.append_child(pugi::node_pcdata).set_value("\n      x y\n    ");
  } else if (space_dimension == 3) {
    names.append_child(pugi::node_pcdata).set_value("\n      x y z\n    ");
  } else if (space_dimension == 4) {
    names.append_child(pugi::node_pcdata).set_value("\n      x y z t\n    ");
  }
}
}  // namespace splinelib::src::io
