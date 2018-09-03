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
  explicit XMLWriterBSpline(std::shared_ptr<spl::BSpline<DIM>> b_spline) {
    this->b_spline = b_spline;
  }

 private:
  double GetDegree(int dimension) override {
    return b_spline->GetDegree(dimension);
  }

  baf::KnotVector GetKnotVector(int dimension) override {
    return b_spline->GetKnotVector(dimension);
  }

  char GetNumberOfControlPoints() override {
    return static_cast<char>(b_spline->GetNumberOfControlPoints());
  }

  char GetSpaceDimension() override {
    return static_cast<char>(b_spline->GetDimension());
  }

  double GetControlPoint(std::array<int, DIM> indices, int dimension) override {
    return b_spline->GetControlPoint(indices, dimension);
  }

  std::array<int, DIM> GetNumberOfPointsInEachDirection() override {
    return b_spline->GetPointsPerDirection();
  }

  std::shared_ptr<spl::BSpline<DIM>> b_spline;
};
}  // namespace io

#endif  // SRC_IO_XML_WRITER_B_SPLINE_H_
