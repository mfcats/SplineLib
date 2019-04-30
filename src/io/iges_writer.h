/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IO_IGES_WRITER_H_
#define SRC_IO_IGES_WRITER_H_

#include <time.h>

#include <any>
#include <fstream>
#include <string>
#include <vector>

#include "b_spline.h"
#include "nurbs.h"
#include "writer.h"

namespace io {
class IGESWriter : public Writer {
 public:
  IGESWriter() = default;

  void WriteFile(const std::vector<std::any> &splines, const char *filename) const override;

 private:
  std::vector<std::string> GetStartSection() const;

  std::vector<std::string> GetGlobalSection(const std::string &filename,
                                            const std::string &delimiter,
                                            const std::string &endDelimiter,
                                            const std::vector<std::any> &splines) const;
  std::vector<std::string> GetGlobalSectionLayout(const std::string *contents) const;

  std::vector<std::string> GetParameterData(const std::string &delimiter, const std::string &endDelimiter,
                                            const std::any &spline, int entityPosition, int *pLine) const;
  std::vector<std::string> GetParameterSectionLayout(const std::string *contents, int entityPosition, int *pLine) const;

  void GetParameterData1D(std::string *contents, const std::string &delimiter, const std::any &spline) const;
  void GetParameterData2D(std::string *contents, const std::string &delimiter, const std::any &spline) const;

  template<int DIM>
  double Get3DControlPoint(std::shared_ptr<spl::Spline<DIM>> spline, int index, int direction) const;

  std::vector<std::string> GetDataEntry(int paramStart, int paramLength, const std::any &spline, int *dLine) const;
  std::vector<std::string> GetDataEntrySectionLayout(const std::string *contents, int *dLine) const;

  std::vector<std::string> GetTerminateSection(int linesS, int linesG, int linesD, int linesP) const;

  int GetDimension(const std::any &spline) const;

  bool IsRational(std::any spline) const;

  void AddToContents(std::string *contents, const std::vector<std::string> &add, const std::string &delimiter) const;

  std::string GetHollerithFormat(const std::string &str) const;

  std::string GetBlock(const std::string &str, int width, bool right) const;

  template<typename T>
  std::string GetString(const T value) const;

  void WriteFile(std::ofstream &file, const std::vector<std::string> &start,
                 const std::vector<std::string> &global, const std::vector<std::string> &data,
                 const std::vector<std::string> &parameter, const std::vector<std::string> &terminate) const;

  void AppendToFile(std::ofstream &file, const std::vector<std::string> &contents) const;

  std::string GetTime() const;

  double GetHighestValue(std::vector<std::any> splines) const;
};
}  // namespace io

#endif  // SRC_IO_IGES_WRITER_H_
