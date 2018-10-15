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

#include "any_casts.h"
#include "b_spline.h"
#include "irit_writer_utils.h"
#include "nurbs.h"

namespace io {
class IRITWriter {
 public:
  IRITWriter() = default;

  void WriteFile(const std::vector<std::any> &splines, const std::string &filename) const {
    std::ofstream newFile;
    newFile.open(filename);
    if (newFile.is_open()) {
      newFile << "[OBJECT SPLINES\n";
      for (unsigned int i = 0; i < splines.size(); i++) {
        AddSpline(newFile, splines[i], i);
        if (i < splines.size() - 1) newFile << "\n";
      }
      newFile << "]";
      newFile.close();
    }
  }

 private:
  void AddSpline(std::ofstream &file, const std::any &spline, int spline_number) const {
    file << "  [OBJECT SPLINE" + std::to_string(spline_number + 1) + "\n    ";
    int spline_dimension = util::AnyCasts::GetSplineDimension(spline);
    switch (spline_dimension) {
      case 1: {
        Write1DSpline(file, spline);
        break;
      }
      case 2: {
        Write2DSpline(file, spline);
        break;
      }
      case 3: {
        Write3DSpline(file, spline);
        break;
      }
      default: {
        throw std::runtime_error("Only splines of dimensions 1 to 3 can be written to an irit file.");
      }
    }
    file << "    ]\n  ]\n";
  }

  void Write1DSpline(std::ofstream &file, const std::any &spline) const {
    std::shared_ptr<spl::Spline<1>> spline_ptr = util::AnyCasts::GetSpline<1>(spline);
    bool rational = util::AnyCasts::IsRational<1>(spline);
    file << "[CURVE BSPLINE "
        + io::IRITWriterUtils<1>::GetNumberOfControlPoints(spline_ptr) + io::IRITWriterUtils<1>::GetOrder(spline_ptr)
        + GetPointType(rational, spline_ptr->GetDimension()) + "\n" + io::IRITWriterUtils<1>::GetKnotVectors(spline_ptr)
        + io::IRITWriterUtils<1>::GetControlPoints(util::AnyCasts::IsRational<1>(spline), spline_ptr, spline);
  }

  void Write2DSpline(std::ofstream &file, const std::any &spline) const {
    std::shared_ptr<spl::Spline<2>> spline_ptr = util::AnyCasts::GetSpline<2>(spline);
    bool rational = util::AnyCasts::IsRational<2>(spline);
    file << "[SURFACE BSPLINE "
        + io::IRITWriterUtils<2>::GetNumberOfControlPoints(spline_ptr) + io::IRITWriterUtils<2>::GetOrder(spline_ptr)
        + GetPointType(rational, spline_ptr->GetDimension()) + "\n" + io::IRITWriterUtils<2>::GetKnotVectors(spline_ptr)
        + io::IRITWriterUtils<2>::GetControlPoints(util::AnyCasts::IsRational<2>(spline), spline_ptr, spline);
  }

  void Write3DSpline(std::ofstream &file, const std::any &spline) const {
    std::shared_ptr<spl::Spline<3>> spline_ptr = util::AnyCasts::GetSpline<3>(spline);
    bool rational = util::AnyCasts::IsRational<3>(spline);
    file << "[TRIVAR BSPLINE "
        + io::IRITWriterUtils<3>::GetNumberOfControlPoints(spline_ptr) + io::IRITWriterUtils<3>::GetOrder(spline_ptr)
        + GetPointType(rational, spline_ptr->GetDimension()) + "\n" + io::IRITWriterUtils<3>::GetKnotVectors(spline_ptr)
        + io::IRITWriterUtils<3>::GetControlPoints(util::AnyCasts::IsRational<3>(spline), spline_ptr, spline);
  }

  std::string GetPointType(bool rational, int space_dimension) const {
    return (rational ? "P" : "E") + std::to_string(space_dimension);
  }
};
}  // namespace io

#endif  // SRC_IO_IRIT_WRITER_H_
