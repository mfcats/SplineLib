/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_IRIT_WRITER_H_
#define SRC_IO_IRIT_WRITER_H_

#include <any>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "b_spline.h"
#include "nurbs.h"

namespace io {
template<int DIM>
class IRITWriter {
 public:
  explicit IRITWriter(const std::vector<std::any> &splines) : splines_(splines) {}

  void WriteIRITFile(const std::string &filename) const {
    std::ofstream newFile;
    newFile.open(filename);
    if (newFile.is_open()) {
      newFile << "[OBJECT SPLINES\n";
      for (uint i = 0; i < splines_.size(); i++) {
        WriteSpline(newFile, i);
      }
      newFile.close();
    } else {
      throw std::runtime_error("The IRIT file couldn't be opened.");  // NOLINT
    }
  }

 private:
  void WriteSpline(std::ofstream &file, uint spline_number) const {
    std::shared_ptr<spl::Spline<DIM>> spline = GetSpline(spline_number);
    file << "  [OBJECT SPLINE" + std::to_string(spline_number + 1) + "\n    [" + GetObjectType() + " BSPLINE "
        + GetNumberOfControlPoints(spline) + GetOrder(spline) + GetPointType(IsRational(spline_number), spline) + "\n"
        + GetKnotVectors(spline) + GetControlPoints(spline_number, IsRational(spline_number), spline) + "    ]\n  ]\n"
        + (spline_number < splines_.size() - 1 ? "\n" : "]");
  }

  std::shared_ptr<spl::Spline<DIM>> GetSpline(int spline_number) const {
    try {
      return std::any_cast<std::shared_ptr<spl::BSpline<DIM>>>(splines_[spline_number]);
    } catch (std::bad_any_cast &msg) {
      try {
        return std::any_cast<std::shared_ptr<spl::NURBS<DIM>>>(splines_[spline_number]);
      } catch (std::bad_any_cast &msg) {
        throw std::runtime_error("Input has to be a vector of pointers to b-splines or nurbs.");
      }
    }
  }

  bool IsRational(int spline_number) const {
    try {
      std::any_cast<std::shared_ptr<spl::NURBS<DIM>>>(splines_[spline_number]);
      return true;
    } catch (std::bad_any_cast &msg) {
      return false;
    }
  }

  std::string GetObjectType() const {
    return DIM == 1 ? "CURVE" : (DIM == 2 ? "SURFACE" : (DIM == 3 ? "TRIVAR" : ""));
  }

  std::string GetNumberOfControlPoints(std::shared_ptr<spl::Spline<DIM>> spline) const {
    std::array<int, DIM> number_of_points = spline->GetPointsPerDirection();
    std::string string;
    for (const auto &points : number_of_points) {
      string += std::to_string(points) + " ";
    }
    return string;
  }

  std::string GetOrder(std::shared_ptr<spl::Spline<DIM>> spline) const {
    std::string string;
    for (int i = 0; i < DIM; i++) {
      string += std::to_string(spline->GetDegree(i).get() + 1) + " ";
    }
    return string;
  }

  std::string GetPointType(bool rational, std::shared_ptr<spl::Spline<DIM>> spline) const {
    return (rational ? "P" : "E") + std::to_string(spline->GetDimension());
  }

  std::string GetKnotVectors(std::shared_ptr<spl::Spline<DIM>> spline) const {
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

  std::string GetControlPoints(int spline_number, bool rational, std::shared_ptr<spl::Spline<DIM>> spline) const {
    std::string string;
    util::MultiIndexHandler<DIM> point_handler(spline->GetPointsPerDirection());
    std::shared_ptr<spl::NURBS<DIM>> nurbs;
    if (rational) {
      nurbs = std::any_cast<std::shared_ptr<spl::NURBS<DIM>>>(splines_[spline_number]);
    }
    for (int control_point = 0; control_point < point_handler.Get1DLength(); ++control_point, point_handler++) {
      string += "      [" + (rational ? std::to_string(nurbs->GetWeight(point_handler.GetIndices())) + " " : "");
      for (int dimension = 0; dimension < spline->GetDimension(); dimension++) {
        string += std::to_string(spline->GetControlPoint(point_handler.GetIndices(), dimension))
            + (dimension < spline->GetDimension() - 1 ? " " : "]\n");
      }
    }
    return string;
  }

  std::vector<std::any> splines_;
};
}  // namespace io

#endif  // SRC_IO_IRIT_WRITER_H_
