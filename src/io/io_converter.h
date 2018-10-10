/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_IO_CONVERTER_H_
#define SRC_IO_IO_CONVERTER_H_

#include <vector>

#include "iges_reader.h"
#include "irit_writer.h"
#include "xml_writer.h"
#include "iges_writer.h"

namespace io {
template<int DIM>
class IOConverter {
 public:
  IOConverter() = default;

  void ConvertIGESFileToXMLFile(const char *input_filename, const char *output_filename) {
    io::IGESReader iges_reader;
    std::vector<std::any> input_splines = iges_reader.ReadIGESFile(input_filename);
    std::vector<std::any> output_splines;
    for (const auto &spline : input_splines) {
      if (GetSplineType(spline) == DIM) {
        output_splines.emplace_back(spline);
      }
    }
    io::XMLWriter<DIM> xml_writer;
    xml_writer.WriteXMLFile(output_splines, output_filename);
  }

 private:
  int GetSplineType(std::any spline) {
    try {
      std::any_cast<std::shared_ptr<spl::BSpline<1>>>(spline);
      return 1;
    } catch (std::bad_any_cast &msg) {
      try {
        std::any_cast<std::shared_ptr<spl::NURBS<1>>>(spline);
        return 1;
      } catch (std::bad_any_cast &msg) {
        try {
          std::any_cast<std::shared_ptr<spl::BSpline<2>>>(spline);
          return 2;
        } catch (std::bad_any_cast &msg) {
          std::any_cast<std::shared_ptr<spl::NURBS<2>>>(spline);
          return 2;
        }
      }
    }
  }
};
}  // namespace io

#endif  // SRC_IO_IO_CONVERTER_H_
