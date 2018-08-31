/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_WRITER_NURBS_H_
#define SRC_IO_WRITER_NURBS_H_

#include <string>

#include "pugixml.hpp"

#include "xml_writer_spline.h"

namespace io {
template<int DIM>
class XMLWriterNURBS : public XMLWriterSpline<DIM> {
 public:
  explicit XMLWriterNURBS(spl::WeightedPhysicalSpace<DIM> physical_space, spl::ParameterSpace<DIM> parameter_space) {
    this->physical_space_ptr = std::make_shared<spl::WeightedPhysicalSpace<DIM>>(physical_space);
    this->parameter_space_ptr = std::make_shared<spl::ParameterSpace<DIM>>(parameter_space);
  }

 private:
  char GetNumberOfControlPoints() override {
    return static_cast<char>(physical_space_ptr->GetNumberOfControlPoints());
  }

  char GetSpaceDimension() override {
    return static_cast<char>(physical_space_ptr->GetDimension());
  }

  void AddWeights(pugi::xml_node *spline) override {
    pugi::xml_node weights = spline->append_child("wght");
    std::string string;
    util::MultiIndexHandler<DIM> weight_handler(GetNumberOfPointsInEachDirection());
    for (int i = 0; i < weight_handler.Get1DLength(); ++i, weight_handler++) {
      auto indices = weight_handler.GetIndices();
      string += "\n      ";
      string += this->GetString(GetWeight(indices)) + "  ";
    }
    weights.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  double GetWeight(std::array<int, DIM> indices) {
    return physical_space_ptr->GetWeight(indices);
  }

  double GetControlPoint(std::array<int, DIM> indices, int dimension) override {
    return physical_space_ptr->GetControlPoint(indices).GetValue(dimension);
  }

  std::array<int, DIM> GetNumberOfPointsInEachDirection() override {
    return physical_space_ptr->GetNumberOfPointsInEachDirection();
  }

  std::shared_ptr<spl::WeightedPhysicalSpace<DIM>> physical_space_ptr;
};
}  // namespace io

#endif  // SRC_IO_XML_WRITER_NURBS_H_
