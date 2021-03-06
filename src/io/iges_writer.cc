/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "iges_writer.h"

#include <iomanip>
#include <limits>
#include <sstream>

#include "any_casts.h"
#include "multi_index_handler.h"
#include "system_operations.h"

void io::IGESWriter::WriteFile(const std::vector<std::any> &splines, const char *filename) const {
  std::ofstream newFile;
  newFile.open(filename);
  if (newFile.is_open()) {
    std::vector<std::string> start = GetStartSection();
    std::vector<std::string> global = GetGlobalSection(filename, ",", ";", splines);
    std::vector<std::string> parameter;
    std::vector<std::string> dataEntry;
    int dLine = 0;
    int pLine = 0;
    int entityPosition = 1;
    int paramStart = 1;
    for (auto &spline : splines) {
      std::vector<std::string> paramTemp = GetParameterData(",", ";", spline, entityPosition, &pLine);
      for (auto &param : paramTemp) {
        parameter.emplace_back(param);
      }
      entityPosition += 2;
      std::vector<std::string>
          dataTemp = GetDataEntry(paramStart, static_cast<int>(paramTemp.size()), spline, &dLine);
      for (auto &data : dataTemp) {
        dataEntry.emplace_back(data);
      }
      paramStart += paramTemp.size();
    }
    std::vector<std::string> terminate = GetTerminateSection(static_cast<int>(start.size()),
                                                             static_cast<int>(global.size()),
                                                             static_cast<int>(dataEntry.size()),
                                                             static_cast<int>(parameter.size()));
    WriteFile(newFile, start, global, dataEntry, parameter, terminate);
  }
  newFile.close();
}

std::vector<std::string> io::IGESWriter::GetStartSection() const {
  return {GetBlock("IGES file generated by SplineLib", 72, false) + "S" + GetBlock(GetString(1), 7, true)};
}

std::vector<std::string> io::IGESWriter::GetGlobalSection(const std::string &filename,
                                                          const std::string &delimiter,
                                                          const std::string &endDelimiter,
                                                          const std::vector<std::any> &splines) const {
  std::string contents;
  AddToContents(&contents, {GetHollerithFormat(delimiter), GetHollerithFormat(endDelimiter),
                            GetHollerithFormat("unknown"), filename, GetHollerithFormat("SplineLib"),
                            GetHollerithFormat("pre_release"),
                            GetString(std::numeric_limits<int>::digits),
                            GetString(std::numeric_limits<float>::max_exponent),
                            GetString(std::numeric_limits<float>::digits),
                            GetString(std::numeric_limits<double>::max_exponent),
                            GetString(std::numeric_limits<double>::digits), GetHollerithFormat("unknown"),
                            GetString(1.0),
                            GetString(2), GetHollerithFormat("MM"), GetString(1), GetString(0.001),
                            GetHollerithFormat(GetTime()), GetString(0.000001), GetString(GetHighestValue(splines)),
                            GetHollerithFormat("unknown"), GetHollerithFormat("unknown"), GetString(11),
                            GetString(0), GetHollerithFormat(GetTime())}, delimiter);
  contents += endDelimiter;
  return GetGlobalSectionLayout(&contents);
}

std::vector<std::string> io::IGESWriter::GetGlobalSectionLayout(const std::string *contents) const {
  int line = 0;
  std::vector<std::string> globalSection;
  for (auto i = 0u; i <= (contents->size() - 1) / 72; ++i) {
    globalSection.emplace_back(GetBlock(contents->substr(i * 72, 72), 72, false)
                                   + 'G' + GetBlock(GetString(++line), 7, true));
  }
  return globalSection;
}

std::vector<std::string> io::IGESWriter::GetParameterData(const std::string &delimiter,
                                                          const std::string &endDelimiter,
                                                          const std::any &spline,
                                                          int entityPosition,
                                                          int *pLine) const {
  std::string contents;
  if (GetDimension(spline) == 126) {
    GetParameterData1D(&contents, delimiter, spline);
  } else if (GetDimension(spline) == 128) {
    GetParameterData2D(&contents, delimiter, spline);
  }
  contents += endDelimiter;
  return GetParameterSectionLayout(&contents, entityPosition, pLine);
}

std::vector<std::string> io::IGESWriter::GetParameterSectionLayout(const std::string *contents,
                                                                   int entityPosition,
                                                                   int *pLine) const {
  std::vector<std::string> parameterData;
  for (auto i = 0u; i <= (contents->size() - 1) / 64; i++) {
    parameterData.emplace_back(GetBlock(contents->substr(i * 64, 64), 64, false)
                                   + ' ' + GetBlock(GetString(entityPosition), 7, true)
                                   + 'P' + GetBlock(GetString(++(*pLine)), 7, true));
  }
  return parameterData;
}

void io::IGESWriter::GetParameterData1D(std::string *contents,
                                        const std::string &delimiter,
                                        const std::any &spline) const {
  std::shared_ptr<spl::Spline<1>> spl;
  int isRational = 0;
  if (IsRational(spline)) {
    spl = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(spline);
  } else {
    spl = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(spline);
    isRational = 1;
  }
  AddToContents(contents,
                {GetString(126),
                 GetString(spl->GetKnotVector(0)->GetNumberOfKnots() - spl->GetDegree(0).get() - 2),
                 GetString(spl->GetDegree(0).get()), GetString(0), GetString(0), GetString(isRational),
                 GetString(0)}, delimiter);

  auto knots = *spl->GetKnotVector(0);
  for (auto &knot : knots) {
    contents->append(GetString(knot.get()) + delimiter);
  }
  for (int i = 0; i < spl->GetNumberOfControlPoints(); ++i) {
    contents->append(GetString(spl->GetWeight({i})) + delimiter);
  }
  for (int i = 0; i < spl->GetNumberOfControlPoints(); ++i) {
    for (int j = 0; j < 3; ++j) {
      contents->append(GetString(Get3DControlPoint(spl, i, j)));
      if (i != spl->GetNumberOfControlPoints() - 1 || j != 2) {
        contents->append(delimiter);
      }
    }
  }
}

void io::IGESWriter::GetParameterData2D(std::string *contents,
                                        const std::string &delimiter,
                                        const std::any &spline) const {
  std::shared_ptr<spl::Spline<2>> spl;
  int isRational = 0;
  if (IsRational(spline)) {
    spl = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(spline);
  } else {
    isRational = 1;
    spl = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(spline);
  }
  AddToContents(contents,
                {GetString(128),
                 GetString(spl->GetKnotVector(0)->GetNumberOfKnots() - spl->GetDegree(0).get() - 2),
                 GetString(spl->GetKnotVector(1)->GetNumberOfKnots() - spl->GetDegree(1).get() - 2),
                 GetString(spl->GetDegree(0).get()), GetString(spl->GetDegree(1).get()),
                 GetString(0), GetString(0), GetString(isRational), GetString(0), GetString(0)},
                delimiter);
  auto knots1 = *spl->GetKnotVector(0);
  auto knots2 = *spl->GetKnotVector(1);
  for (auto &i : knots1) {
    contents->append(GetString(i.get()) + delimiter);
  }
  for (auto &i : knots2) {
    contents->append(GetString(i.get()) + delimiter);
  }
  for (int i = 0; i < spl->GetNumberOfControlPoints(); ++i) {
    contents->append(GetString(spl->GetWeight({i})) + delimiter);
  }
  for (int i = 0; i < spl->GetNumberOfControlPoints(); ++i) {
    for (int j = 0; j < 3; ++j) {
      contents->append(GetString(Get3DControlPoint(spl, i, j)));
      if (i != spl->GetNumberOfControlPoints() - 1 || j != 2) {
        contents->append(delimiter);
      }
    }
  }
}

template<int DIM>
double io::IGESWriter::Get3DControlPoint(std::shared_ptr<spl::Spline<DIM>> spline, int index, int direction) const {
  if (direction < spline->GetPointDim()) {
    return spline->GetControlPoint({index}, direction);
  }
  return 0.0;
}

std::vector<std::string> io::IGESWriter::GetDataEntry(int paramStart,
                                                      int paramLength,
                                                      const std::any &spline,
                                                      int *dLine) const {
  std::string contents;
  AddToContents(&contents,
                {GetBlock(GetString(GetDimension(spline)), 8, true), GetBlock(GetString(paramStart), 8, true),
                 GetBlock(GetString(0), 8, true), GetBlock(GetString(0), 8, true), GetBlock(GetString(0), 8, true),
                 GetBlock(GetString(0), 8, true), GetBlock(GetString(0), 8, true), GetBlock(GetString(0), 8, true),
                 " 0 0 0 1", GetBlock(GetString(GetDimension(spline)), 8, true),
                 GetBlock(GetString(0), 8, true), GetBlock(GetString(0), 8, true),
                 GetBlock(GetString(paramLength), 8, true), GetBlock(GetString(0), 8, true),
                 GetBlock("", 16, true), GetBlock("SPLINE", 8, true), GetBlock(GetString(0), 8, true)}, "");
  return GetDataEntrySectionLayout(&contents, dLine);
}

std::vector<std::string> io::IGESWriter::GetDataEntrySectionLayout(const std::string *contents, int *dLine) const {
  std::vector<std::string> dataEntrySection;
  for (auto i = 0u; i <= (contents->size() - 1) / 72; i++) {
    dataEntrySection.emplace_back(GetBlock(contents->substr(i * 72, 72), 72, false)
                                      + 'D' + GetBlock(GetString(++(*dLine)), 7, true));
  }
  return dataEntrySection;
}

std::vector<std::string> io::IGESWriter::GetTerminateSection(int linesS, int linesG, int linesD, int linesP) const {
  return {'S' + GetBlock(GetString(linesS), 7, true) + 'G' + GetBlock(GetString(linesG), 7, true) + 'D'
              + GetBlock(GetString(linesD), 7, true) + 'P' + GetBlock(GetString(linesP), 7, true)
              + GetBlock("T", 41, true) + GetBlock(GetString(1), 7, true)};
}

int io::IGESWriter::GetDimension(const std::any &spline) const {
  if (util::AnyCasts::GetSplineDimension(spline) == 1) {
    return 126;
  }
  if (util::AnyCasts::GetSplineDimension(spline) == 2) {
    return 128;
  }
  throw std::runtime_error("Only splines of dimensions 1 or 2 can be written to an iges file.");
}

bool io::IGESWriter::IsRational(std::any spline) const {
  try {
    std::any_cast<std::shared_ptr<spl::NURBS<1>>>(spline);
    return true;
  } catch (std::bad_any_cast &msg) {
    try {
      std::any_cast<std::shared_ptr<spl::NURBS<2>>>(spline);
      return true;
    } catch (std::bad_any_cast &msg) {
      return false;
    }
  }
}

void io::IGESWriter::AddToContents(std::string *contents,
                                   const std::vector<std::string> &add,
                                   const std::string &delimiter) const {
  for (const auto &i : add) {
    contents->append(i + delimiter);
  }
}

std::string io::IGESWriter::GetHollerithFormat(const std::string &str) const {
  return GetString(str.size()) + "H" + str;
}

std::string io::IGESWriter::GetBlock(const std::string &str, int width, bool right) const {
  return right ? (std::string(width - str.size(), ' ') + str) : (str + std::string(width - str.size(), ' '));
}

template<typename T>
std::string io::IGESWriter::GetString(const T value) const {
  std::ostringstream out;
  out.precision(18);
  out << value;
  return out.str();
}

void io::IGESWriter::WriteFile(std::ofstream &file,
                               const std::vector<std::string> &start,
                               const std::vector<std::string> &global,
                               const std::vector<std::string> &data,
                               const std::vector<std::string> &parameter,
                               const std::vector<std::string> &terminate) const {
  AppendToFile(file, start);
  AppendToFile(file, global);
  AppendToFile(file, data);
  AppendToFile(file, parameter);
  AppendToFile(file, terminate);
}

void io::IGESWriter::AppendToFile(std::ofstream &file, const std::vector<std::string> &contents) const {
  for (auto &content : contents) {
    file << content << std::endl;
  }
}

std::string io::IGESWriter::GetTime() const {
  struct tm ptm = util::SystemOperations::GetTime();
  std::string date = GetString((ptm.tm_year + 1900) * 10000 + (ptm.tm_mon + 1) * 100 + ptm.tm_mday);
  std::string time = GetString(ptm.tm_hour * 10000 + ptm.tm_min * 100 + ptm.tm_sec);
  return GetHollerithFormat(date + "." + time);
}

double io::IGESWriter::GetHighestValue(std::vector<std::any> splines) const {
  double highestValue = 0;
  for (auto &spline : splines) {
    int dimension = GetDimension(spline);
    bool isRational = IsRational(spline);
    if (dimension == 126) {
      if (isRational) {
        auto spl = std::any_cast<std::shared_ptr<spl::NURBS<1>>>(spline);
        if (highestValue < spl->GetExpansion()) {
          highestValue = spl->GetExpansion();
        }
      } else {
        auto spl = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(spline);
        if (highestValue < spl->GetExpansion()) {
          highestValue = spl->GetExpansion();
        }
      }
    } else if (dimension == 128) {
      if (isRational) {
        auto spl = std::any_cast<std::shared_ptr<spl::NURBS<2>>>(spline);
        if (highestValue < spl->GetExpansion()) {
          highestValue = spl->GetExpansion();
        }
      } else {
        auto spl = std::any_cast<std::shared_ptr<spl::BSpline<2>>>(spline);
        if (highestValue < spl->GetExpansion()) {
          highestValue = spl->GetExpansion();
        }
      }
    }
  }
  return highestValue;
}
