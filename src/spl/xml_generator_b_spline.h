/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_XML_GENERATOR_B_SPLINE_H_
#define SRC_SPL_XML_GENERATOR_B_SPLINE_H_

#include "pugixml.hpp"

#include "physical_space.h"
#include "xml_generator_spline.h"

namespace spl {
template<int DIM>
class XMLGenerator_B_Spline : public XMLGenerator_Spline<DIM> {
 public:
  explicit XMLGenerator_B_Spline(PhysicalSpace<DIM> physical_space, ParameterSpace<DIM> parameter_space) {
    this->physical_space_ptr = std::make_shared<PhysicalSpace<DIM>>(physical_space);
    this->parameter_space_ptr = std::make_shared<ParameterSpace<DIM>>(parameter_space);
  }

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

 private:
  std::shared_ptr<PhysicalSpace<DIM>> physical_space_ptr;
};
}  // namespace spl

#endif  // SRC_SPL_XML_GENERATOR_B_SPLINE_H_
