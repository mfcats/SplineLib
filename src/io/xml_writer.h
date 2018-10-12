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
#include <string>
#include <vector>

#include "pugixml.hpp"

#include "any_casts.h"
#include "b_spline.h"
#include "nurbs.h"
#include "xml_writer_utils.h"

namespace io {
class XMLWriter {
 public:
  XMLWriter() = default;

  void WriteXMLFile(const std::vector<std::any> &splines, const char *filename) {
    pugi::xml_document doc;
    pugi::xml_node spline_list = doc.append_child("SplineList");
    spline_list.append_attribute("NumberOfSplines") = std::to_string(splines.size()).c_str();
    for (const auto &spline : splines) {
      AddSpline(&spline_list, spline);
    }
    doc.save_file(filename, "  ", pugi::format_indent_attributes, pugi::encoding_utf8);
  }

 private:
  int GetSplineDimension(const std::any &spline) {
    try {
      std::any_cast<std::shared_ptr<spl::BSpline<1>>>(spline);
      return 1;
    } catch (std::bad_any_cast &msg) {
      try {
        std::any_cast<std::shared_ptr<spl::NURBS<1>>>(spline);
        return 1;
      } catch (std::bad_any_cast &msg) {
        try {
          std::any_cast<std::shared_ptr<spl::BSpline<2>>>(spline);
          return 2;
        } catch (std::bad_any_cast &msg) {
          try {
            std::any_cast<std::shared_ptr<spl::NURBS<2>>>(spline);
            return 2;
          } catch (std::bad_any_cast &msg) {
            try {
              std::any_cast<std::shared_ptr<spl::BSpline<3>>>(spline);
              return 3;
            } catch (std::bad_any_cast &msg) {
              try {
                std::any_cast<std::shared_ptr<spl::NURBS<3>>>(spline);
                return 3;
              } catch (std::bad_any_cast &msg) {
                try {
                  std::any_cast<std::shared_ptr<spl::BSpline<4>>>(spline);
                  return 4;
                } catch (std::bad_any_cast &msg) {
                  try {
                    std::any_cast<std::shared_ptr<spl::NURBS<4>>>(spline);
                    return 4;
                  } catch (std::bad_any_cast &msg) {

                    throw std::runtime_error(
                        "Input has to be a pointer to a b-spline or nurbs of dimension 1, 2, 3 or 4.");
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  void AddSpline(pugi::xml_node *spline_list, const std::any &spline) {
    int spline_dimension = GetSplineDimension(spline);
    switch (spline_dimension) {
      case 1: {
        Add1DSpline(spline_list, spline);
        break;
      }
      case 2: {
        Add2DSpline(spline_list, spline);
        break;
      }
      case 3: {
        Add3DSpline(spline_list, spline);
        break;
      }
      case 4: {
        Add4DSpline(spline_list, spline);
        break;
      }
      default: {
        throw std::runtime_error("Only splines of dimensions 1 to 4 can be written into a xml file.");
      }
    }
  }

  void Add1DSpline(pugi::xml_node *spline_list, const std::any &spline) {
    std::shared_ptr<spl::Spline<1>> spline_ptr = util::AnyCasts<1>::GetSpline(spline);
    pugi::xml_node spline_node = spline_list->append_child("SplineEntry");
    AddSplineAttributes(&spline_node, 1, spline_ptr->GetDimension(), spline_ptr->GetNumberOfControlPoints());
    AddControlPointVarNames(&spline_node, spline_ptr->GetDimension());
    io::XMLWriterUtils<1>::AddControlPointVars(&spline_node, spline_ptr);
    if (util::AnyCasts<1>::IsRational(spline)) io::XMLWriterUtils<1>::AddWeights(&spline_node, spline);
    io::XMLWriterUtils<1>::AddDegrees(&spline_node, spline_ptr);
    io::XMLWriterUtils<1>::AddKnotVectors(&spline_node, spline_ptr);
  }

  void Add2DSpline(pugi::xml_node *spline_list, const std::any &spline) {
    std::shared_ptr<spl::Spline<2>> spline_ptr = util::AnyCasts<2>::GetSpline(spline);
    pugi::xml_node spline_node = spline_list->append_child("SplineEntry");
    AddSplineAttributes(&spline_node, 2, spline_ptr->GetDimension(), spline_ptr->GetNumberOfControlPoints());
    AddControlPointVarNames(&spline_node, spline_ptr->GetDimension());
    io::XMLWriterUtils<2>::AddControlPointVars(&spline_node, spline_ptr);
    if (util::AnyCasts<2>::IsRational(spline)) io::XMLWriterUtils<2>::AddWeights(&spline_node, spline);
    io::XMLWriterUtils<2>::AddDegrees(&spline_node, spline_ptr);
    io::XMLWriterUtils<2>::AddKnotVectors(&spline_node, spline_ptr);
  }

  void Add3DSpline(pugi::xml_node *spline_list, const std::any &spline) {
    std::shared_ptr<spl::Spline<3>> spline_ptr = util::AnyCasts<3>::GetSpline(spline);
    pugi::xml_node spline_node = spline_list->append_child("SplineEntry");
    AddSplineAttributes(&spline_node, 3, spline_ptr->GetDimension(), spline_ptr->GetNumberOfControlPoints());
    AddControlPointVarNames(&spline_node, spline_ptr->GetDimension());
    io::XMLWriterUtils<3>::AddControlPointVars(&spline_node, spline_ptr);
    if (util::AnyCasts<3>::IsRational(spline)) io::XMLWriterUtils<3>::AddWeights(&spline_node, spline);
    io::XMLWriterUtils<3>::AddDegrees(&spline_node, spline_ptr);
    io::XMLWriterUtils<3>::AddKnotVectors(&spline_node, spline_ptr);
  }

  void Add4DSpline(pugi::xml_node *spline_list, const std::any &spline) {
    std::shared_ptr<spl::Spline<4>> spline_ptr = util::AnyCasts<4>::GetSpline(spline);
    pugi::xml_node spline_node = spline_list->append_child("SplineEntry");
    AddSplineAttributes(&spline_node, 4, spline_ptr->GetDimension(), spline_ptr->GetNumberOfControlPoints());
    AddControlPointVarNames(&spline_node, spline_ptr->GetDimension());
    io::XMLWriterUtils<4>::AddControlPointVars(&spline_node, spline_ptr);
    if (util::AnyCasts<4>::IsRational(spline)) io::XMLWriterUtils<4>::AddWeights(&spline_node, spline);
    io::XMLWriterUtils<4>::AddDegrees(&spline_node, spline_ptr);
    io::XMLWriterUtils<4>::AddKnotVectors(&spline_node, spline_ptr);
  }

  void AddSplineAttributes(pugi::xml_node *spline_node, int spline_dimension, int space_dimension, int control_points) {
    spline_node->append_attribute("splDim") = spline_dimension;
    spline_node->append_attribute("spaceDim") = space_dimension;
    spline_node->append_attribute("numOfCntrlPntVars") = space_dimension;
    spline_node->append_attribute("numCntrlPnts") = control_points;
  }

  void AddControlPointVarNames(pugi::xml_node *spline_node, int space_dimension) {
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
};
}  // namespace io

#endif  // SRC_IO_XML_WRITER_H_
