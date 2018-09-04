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
  explicit XMLWriterBSpline(std::vector<spl::BSpline<DIM>> b_splines)
      : XMLWriterSpline<DIM>(static_cast<int>(b_splines.size())) {
    for (auto &b_spline : b_splines) {
      this->b_splines.push_back(std::make_shared<spl::BSpline<DIM>>(b_spline));
    }
  }

 private:
  double GetDegree(int spline, int dimension) override {
    return b_splines[spline]->GetDegree(dimension);
  }

  baf::KnotVector GetKnotVector(int spline, int dimension) override {
    return b_splines[spline]->GetKnotVector(dimension);
  }

  char GetNumberOfControlPoints(int spline) override {
    return static_cast<char>(b_splines[spline]->GetNumberOfControlPoints());
  }

  char GetSpaceDimension(int spline) override {
    return static_cast<char>(b_splines[spline]->GetDimension());
  }

  double GetControlPoint(int spline, std::array<int, DIM> indices, int dimension) override {
    return b_splines[spline]->GetControlPoint(indices, dimension);
  }

  std::array<int, DIM> GetNumberOfPointsInEachDirection(int spline) override {
    return b_splines[spline]->GetPointsPerDirection();
  }

  std::vector<std::shared_ptr<spl::BSpline<DIM>>> b_splines;
};
}  // namespace io

#endif  // SRC_IO_XML_WRITER_B_SPLINE_H_
