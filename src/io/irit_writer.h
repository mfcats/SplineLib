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

namespace io {
template<int DIM>
class IRITWriter {
 public:
  explicit IRITWriter(std::vector<std::any> splines) {
    for (auto &spline : splines) {
      this->splines_.push_back(std::make_shared<std::any>(spline));
    }
  }

  void WriteIRITFile(const std::string &filename) {
    std::ofstream newFile;
    newFile.open(filename);
    if (newFile.is_open()) {
      std::string object_type = GetObjectType();
      AppendFile(newFile, {"[", object_type, " BSPLINE "});
      std::vector<std::string>
          number_of_control_points = GetNumberOfControlPoints(std::any_cast<spl::BSpline<DIM>>(*splines_[0]));
      AppendFile(newFile, number_of_control_points);
      std::vector<std::string> order = GetOrder(std::any_cast<spl::BSpline<DIM>>(*splines_[0]));
      AppendFile(newFile, order);
      std::string point_type = GetPointType(false, std::any_cast<spl::BSpline<DIM>>(*splines_[0]));
      AppendFile(newFile, {point_type, "\n"});
      std::vector<std::string> knot_vectors = GetKnotVectors(std::any_cast<spl::BSpline<DIM>>(*splines_[0]));
      AppendFile(newFile, knot_vectors);
      std::vector<std::string> control_points = GetControlPoints(std::any_cast<spl::BSpline<DIM>>(*splines_[0]));
      AppendFile(newFile, control_points);
      AppendFile(newFile, {"]"});
      newFile.close();
    } else {
      throw std::runtime_error("The IRIT file couldn't be opened.");  // NOLINT
    }
  }

  void AppendFile(std::ofstream &file, const std::vector<std::string> &contents) {
    for (auto &content : contents) {
      file << content;
    }
  }

 private:
  std::string GetObjectType() {
    return DIM == 1 ? "CURVE" : (DIM == 2 ? "SURFACE" : (DIM == 3 ? "TRIVAR" : ""));
  }

  std::vector<std::string> GetNumberOfControlPoints(spl::BSpline<DIM> spline) {
    std::array<int, DIM> number = spline.GetPointsPerDirection();
    std::vector<std::string> string_vector;
    for (const auto &points : number) {
      string_vector.emplace_back(std::to_string(points));
      string_vector.emplace_back(" ");
    }
    return string_vector;
  }

  std::vector<std::string> GetOrder(spl::BSpline<DIM> spline) {
    std::vector<std::string> string_vector;
    for (int i = 0; i < DIM; i++) {
      string_vector.emplace_back(std::to_string(spline.GetDegree(i).get() + 1));
      string_vector.emplace_back(" ");
    }
    return string_vector;
  }

  std::string GetPointType(bool rational, spl::BSpline<DIM> spline) {
    return (rational ? "P" : "E") + std::to_string(spline.GetDimension());
  }

  std::vector<std::string> GetKnotVectors(spl::BSpline<DIM> spline) {
    std::vector<std::string> string_vector;
    for (int i = 0; i < DIM; i++) {
      string_vector.emplace_back("[KV ");
      std::shared_ptr<baf::KnotVector> knot_vector = spline.GetKnotVector(i);
      for (const auto &knot : *knot_vector) {
        string_vector.emplace_back(std::to_string(knot.get()));
        string_vector.emplace_back(" ");
      }
      string_vector.emplace_back("]\n");
    }
    return string_vector;
  }

  std::vector<std::string> GetControlPoints(spl::BSpline<DIM> spline) {
    std::vector<std::string> string_vector;
    util::MultiIndexHandler<DIM> control_point_handler(spline.GetPointsPerDirection());
    for (int i = 0; i < control_point_handler.Get1DLength(); ++i, control_point_handler++) {
      string_vector.push_back("[");
      for (int j = 0; j < spline.GetDimension(); j++) {
        string_vector.push_back(std::to_string(spline.GetControlPoint(control_point_handler.GetIndices(), j)) + " ");
      }
      string_vector.push_back("]\n");
    }
    return string_vector;
  }

  std::vector<std::shared_ptr<std::any>> splines_;
};
}  // namespace io

#endif  // SRC_IO_IRIT_WRITER_H_
