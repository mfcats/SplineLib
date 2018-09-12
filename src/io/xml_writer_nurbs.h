/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_XML_WRITER_NURBS_H_
#define SRC_IO_XML_WRITER_NURBS_H_

#include <string>
#include <vector>

#include "pugixml.hpp"

#include "xml_writer_spline.h"

namespace io {
template<int DIM>
class XMLWriterNURBS : public XMLWriterSpline<DIM> {
 public:
  explicit XMLWriterNURBS(std::vector<spl::NURBS<DIM>> nurbs)
      : XMLWriterSpline<DIM>(static_cast<int>(nurbs.size())) {
    for (auto &spline : nurbs) {
      this->nurbs.push_back(std::make_shared<spl::NURBS<DIM>>(spline));
    }
  }

 private:
  Degree GetDegree(int spline, int dimension) override {
    return nurbs[spline]->GetDegree(dimension);
  }

  std::shared_ptr<baf::KnotVector> GetKnotVector(int spline, int dimension) override {
    return nurbs[spline]->GetKnotVector(dimension);
  }

  char GetNumberOfControlPoints(int spline) override {
    return static_cast<char>(nurbs[spline]->GetNumberOfControlPoints());
  }

  char GetSpaceDimension(int spline) override {
    return static_cast<char>(nurbs[spline]->GetDimension());
  }

  void AddWeights(pugi::xml_node *spline, int number) override {
    pugi::xml_node weights = spline->append_child("wght");
    std::string string;
    util::MultiIndexHandler<DIM> weight_handler(GetNumberOfPointsInEachDirection(number));
    for (int i = 0; i < weight_handler.Get1DLength(); ++i, weight_handler++) {
      auto indices = weight_handler.GetIndices();
      string += "\n      " + std::to_string(GetWeight(number, indices)) + "  ";
    }
    weights.append_child(pugi::node_pcdata).text() = (string + "\n    ").c_str();
  }

  double GetWeight(int spline, std::array<int, DIM> indices) {
    return nurbs[spline]->GetWeight(indices);
  }

  double GetControlPoint(int spline, std::array<int, DIM> indices, int dimension) override {
    return nurbs[spline]->GetControlPoint(indices, dimension);
  }

  std::array<int, DIM> GetNumberOfPointsInEachDirection(int spline) override {
    return nurbs[spline]->GetPointsPerDirection();
  }

  std::vector<std::shared_ptr<spl::NURBS<DIM>>> nurbs;
};
}  // namespace io

#endif  // SRC_IO_XML_WRITER_NURBS_H_
