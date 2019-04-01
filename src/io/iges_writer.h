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

namespace io {
class IGESWriter {
 public:
  IGESWriter() = default;

  void WriteFile(std::vector<std::any> splines, const std::string &filename);

 private:
  std::vector<std::string> GetStartSection();

  std::vector<std::string> GetGlobalSection(const std::string &filename, const std::string &delimiter,
                                            const std::string &endDelimiter, const std::vector<std::any> &splines);
  std::vector<std::string> GetGlobalSectionLayout(const std::string *contents);

  std::vector<std::string> GetParameterData(const std::string &delimiter, const std::string &endDelimiter,
                                            const std::any &spline, int entityPosition, int *pLine);
  std::vector<std::string> GetParameterSectionLayout(const std::string *contents, int entityPosition, int *pLine);

  void GetParameterData1D(std::string *contents, const std::string &delimiter, const std::any &spline);
  void GetParameterData2D(std::string *contents, const std::string &delimiter, const std::any &spline);

  template<int DIM>
  double Get3DControlPoint(std::shared_ptr<spl::Spline<DIM>> spline, int index, int direction);

  std::vector<std::string> GetDataEntry(int paramStart, int paramLength, const std::any &spline, int *dLine);
  std::vector<std::string> GetDataEntrySectionLayout(const std::string *contents, int *dLine);

  std::vector<std::string> GetTerminateSection(int linesS, int linesG, int linesD, int linesP);

  int GetDimension(const std::any &spline);

  bool IsRational(std::any spline);

  void AddToContents(std::string *contents, const std::vector<std::string> &add, const std::string &delimiter);

  std::string GetHollerithFormat(const std::string &str);

  std::string GetBlock(const std::string &str, int width, bool right);

  template<typename T>
  std::string GetString(const T value);

  void WriteFile(std::ofstream &file, const std::vector<std::string> &start,
                 const std::vector<std::string> &global, const std::vector<std::string> &data,
                 const std::vector<std::string> &parameter, const std::vector<std::string> &terminate);

  void AppendToFile(std::ofstream &file, const std::vector<std::string> &contents);

  std::string GetTime();

  double GetHighestValue(std::vector<std::any> splines);
};
}  // namespace io

#endif  // SRC_IO_IGES_WRITER_H_
