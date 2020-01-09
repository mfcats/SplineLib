/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_IO_IRIT_READER_H_
#define SRC_IO_IRIT_READER_H_

#include <any>
#include <string>
#include <vector>

#include "src/spl/b_spline.h"
#include "reader.h"

namespace splinelib::src::io {
class IRITReader : public Reader {
 public:
  IRITReader() = default;

  std::vector<std::any> ReadFile(const char *filename) override;

 private:
  std::vector<int> GetSplinePositions(const std::vector<std::string> &entries) const;

  static int GetDimension(const std::string &type);

  std::any Get1DSpline(int start, const std::vector<std::string> &entries) const;
  std::any Get2DSpline(int start, const std::vector<std::string> &entries) const;
  std::any Get3DSpline(int start, const std::vector<std::string> &entries) const;

  static int GetTotalNumberOfControlPoints(int start, const std::vector<std::string> &entries);

  std::vector<spl::ControlPoint> GetControlPoints(int start, const std::vector<std::string> &entries,
                                                  bool rational, std::vector<Weight> weights) const;

  std::vector<Weight> GetWeights(int start, const std::vector<std::string> &entries, bool rational) const;

  int GetPositionOfFirstControlPoint(int start, const std::vector<std::string> &entries) const;
};
}  // namespace splinelib::src::io

#endif  // SRC_IO_IRIT_READER_H_
