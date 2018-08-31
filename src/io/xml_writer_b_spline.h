/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_XML_WRITER_B_SPLINE_H_
#define SRC_IO_XML_WRITER_B_SPLINE_H_

#include "pugixml.hpp"

#include "xml_writer_spline.h"

namespace io {
template<int DIM>
class XMLWriterBSpline : public XMLWriterSpline<DIM> {
 public:
  XMLWriterBSpline(spl::PhysicalSpace<DIM> physical_space, spl::ParameterSpace<DIM> parameter_space) {
    this->physical_space_ptr = std::make_shared<spl::PhysicalSpace<DIM>>(physical_space);
    this->parameter_space_ptr = std::make_shared<spl::ParameterSpace<DIM>>(parameter_space);
  }

 private:
  char GetNumberOfControlPoints() override {
    return static_cast<char>(physical_space_ptr->GetNumberOfControlPoints());
  }

  char GetSpaceDimension() override {
    return static_cast<char>(physical_space_ptr->GetDimension());
  }

  double GetControlPoint(std::array<int, DIM> indices, int dimension) override {
    return physical_space_ptr->GetControlPoint(indices).GetValue(dimension);
  }

  std::array<int, DIM> GetNumberOfPointsInEachDirection() override {
    return physical_space_ptr->GetNumberOfPointsInEachDirection();
  }

  std::shared_ptr<spl::PhysicalSpace<DIM>> physical_space_ptr;
};
}  // namespace io

#endif  // SRC_IO_XML_WRITER_B_SPLINE_H_
