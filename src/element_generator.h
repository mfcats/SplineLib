/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_ELEMENT_GENERATOR_H
#define SPLINELIB_ELEMENT_GENERATOR_H

#include "b_spline.h"
#include "element.h"

class ElementGenerator {
 public:
  ElementGenerator(int degree_, KnotVector knot_vector_, std::vector<ControlPoint> control_points_);

  std::vector<Element> GetElementList() {
    std::vector<Element> elements;
    for (uint64_t currentKnotSpan = 0;
         currentKnotSpan < knot_vector_.Size() - degree_ - 1; currentKnotSpan++) {
      double upper = knot_vector_.knot(currentKnotSpan + 1);
      double lower = knot_vector_.knot(currentKnotSpan);
      if ((upper - lower) != 0) {
        elements.emplace_back(Element(1, control_points_));
      }
    }
  }

 private:
  int degree_;
  KnotVector knot_vector_;
  std::vector<ControlPoint> control_points_;
};

#endif //SPLINELIB_ELEMENT_GENERATOR_H
