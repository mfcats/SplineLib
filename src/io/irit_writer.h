/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_IO_IRIT_WRITER_H_
#define SRC_IO_IRIT_WRITER_H_

#include <any>
#include <fstream>
#include <string>
#include <vector>

#include "writer.h"

namespace splinelib::src::io {
class IRITWriter : public Writer {
 public:
  IRITWriter() = default;

  void WriteFile(const std::vector<std::any> &splines, const char *filename) const override;

 private:
  void AddSpline(std::ofstream &file, const std::any &spline, int spline_number) const;

  void Write1DSpline(std::ofstream &file, const std::any &spline) const;
  void Write2DSpline(std::ofstream &file, const std::any &spline) const;
  void Write3DSpline(std::ofstream &file, const std::any &spline) const;

  std::string GetPointType(bool rational, int space_dimension) const;
};
}  // namespace splinelib::src::io

#endif  // SRC_IO_IRIT_WRITER_H_
