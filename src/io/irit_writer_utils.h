/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_IRIT_WRITER_UTILS_H_
#define SRC_IO_IRIT_WRITER_UTILS_H_

#include <string>

#include "any_casts.h"
#include "b_spline.h"
#include "nurbs.h"

namespace io {
template<int DIM>
class IRITWriterUtils {
 public:
  static std::string GetNumberOfControlPoints(std::shared_ptr<spl::Spline<DIM>> spline) {
    std::array<int, DIM> number_of_points = spline->GetPointsPerDirection();
    std::string string;
    for (const auto &points : number_of_points) {
      string += std::to_string(points) + " ";
    }
    return string;
  }

  static std::string GetOrder(std::shared_ptr<spl::Spline<DIM>> spline) {
    std::string string;
    for (int i = 0; i < DIM; i++) {
      string += std::to_string(spline->GetDegree(i).get() + 1) + " ";
    }
    return string;
  }

  static std::string GetKnotVectors(std::shared_ptr<spl::Spline<DIM>> spline) {
    std::string string;
    for (int dimension = 0; dimension < DIM; dimension++) {
      string += "      [KV ";
      std::shared_ptr<baf::KnotVector> knot_vector = spline->GetKnotVector(dimension);
      for (size_t knot = 0; knot < knot_vector->GetNumberOfKnots(); knot++) {
        string += std::to_string(knot_vector->GetKnot(knot).get()) +
            (knot < knot_vector->GetNumberOfKnots() - 1 ? " " : "]\n");
      }
    }
    return string;
  }

  static std::string
  GetControlPoints(bool rational, std::shared_ptr<spl::Spline<DIM>> spline_ptr, const std::any &spline) {
    std::string string;
    util::MultiIndexHandler<DIM> point_handler(spline_ptr->GetPointsPerDirection());
    std::shared_ptr<spl::NURBS<DIM>> nurbs;
    if (rational) {
      nurbs = std::any_cast<std::shared_ptr<spl::NURBS<DIM>>>(spline);
    }
    for (int control_point = 0; control_point < point_handler.Get1DLength(); ++control_point, point_handler++) {
      string += "      [" + (rational ? std::to_string(nurbs->GetWeight(point_handler.GetIndices())) + " " : "");
      for (int dimension = 0; dimension < spline_ptr->GetDimension(); dimension++) {
        string += std::to_string(spline_ptr->GetControlPoint(point_handler.GetIndices(), dimension)
                                     * (rational ? nurbs->GetWeight(point_handler.GetIndices()) : 1))
            + (dimension < spline_ptr->GetDimension() - 1 ? " " : "]\n");
      }
    }
    return string;
  }
};
}  // namespace io

#endif  // SRC_IO_IRIT_WRITER_UTILS_H_
