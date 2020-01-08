/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "src/io/iges_writer.h"

#include <iomanip>
#include <limits>
#include <sstream>

#include "src/util/any_casts.h"
#include "src/util/multi_index_handler.h"
#include "src/util/string_operations.h"
#include "src/util/system_operations.h"

namespace splinelib::src::io {
void IGESWriter::WriteFile(const std::vector<std::any> &splines, const char *filename) const {
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

std::vector<std::string> IGESWriter::GetStartSection() const {
  return {GetBlock("IGES file generated by SplineLib", 72, false) + "S" + GetBlock("1", 7, true)};
}

std::vector<std::string> IGESWriter::GetGlobalSection(const std::string &filename,
                                                          const std::string &delimiter,
                                                          const std::string &endDelimiter,
                                                          const std::vector<std::any> &splines) const {
  std::string contents;
  AddToContents(&contents, {GetHollerithFormat(delimiter), GetHollerithFormat(endDelimiter),
                            GetHollerithFormat("unknown"), filename, GetHollerithFormat("SplineLib"),
                            GetHollerithFormat("pre_release"), std::to_string(std::numeric_limits<int>::digits),
                            std::to_string(std::numeric_limits<float>::max_exponent),
                            std::to_string(std::numeric_limits<float>::digits),
                            std::to_string(std::numeric_limits<double>::max_exponent),
                            std::to_string(std::numeric_limits<double>::digits), GetHollerithFormat("unknown"), "1.0",
                            "2", GetHollerithFormat("MM"), "1", "0.001", GetHollerithFormat(GetTime()), "0.000001",
                            util::string_operations::GetStringWithHighPrecision(GetHighestValue(splines)),
                            GetHollerithFormat("unknown"), GetHollerithFormat("unknown"), "11", "0",
                            GetHollerithFormat(GetTime())}, delimiter);
  contents += endDelimiter;
  return GetGlobalSectionLayout(&contents);
}

std::vector<std::string> IGESWriter::GetGlobalSectionLayout(const std::string *contents) const {
  int line = 0;
  std::vector<std::string> globalSection;
  for (auto i = 0u; i <= (contents->size() - 1) / 72; ++i) {
    globalSection.emplace_back(GetBlock(contents->substr(i * 72, 72), 72, false)
                                   + 'G' + GetBlock(std::to_string(++line), 7, true));
  }
  return globalSection;
}

std::vector<std::string> IGESWriter::GetParameterData(const std::string &delimiter,
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

std::vector<std::string> IGESWriter::GetParameterSectionLayout(const std::string *contents,
                                                                   int entityPosition,
                                                                   int *pLine) const {
  std::vector<std::string> parameterData;
  for (auto i = 0u; i <= (contents->size() - 1) / 64; i++) {
    parameterData.emplace_back(GetBlock(contents->substr(i * 64, 64), 64, false)
                                   + ' ' + GetBlock(std::to_string(entityPosition), 7, true)
                                   + 'P' + GetBlock(std::to_string(++(*pLine)), 7, true));
  }
  return parameterData;
}

void IGESWriter::GetParameterData1D(std::string *contents,
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
                {"126", std::to_string(spl->GetKnotVector(0)->GetNumberOfKnots() - spl->GetDegree(0).Get() - 2),
                 std::to_string(spl->GetDegree(0).Get()), "0", "0", std::to_string(isRational), "0"}, delimiter);

  auto knots = *spl->GetKnotVector(0);
  for (auto &knot : knots) {
    contents->append(util::string_operations::GetStringWithHighPrecision(knot.Get()) + delimiter);
  }
  for (int i = 0; i < spl->GetTotalNumberOfControlPoints(); ++i) {
    contents->append(util::string_operations::GetStringWithHighPrecision(spl->GetWeight(i)) + delimiter);
  }
  for (int i = 0; i < spl->GetTotalNumberOfControlPoints(); ++i) {
    for (int j = 0; j < 3; ++j) {
      contents->append(util::string_operations::GetStringWithHighPrecision(Get3DControlPoint(spl, i, j)));
      if (i != spl->GetTotalNumberOfControlPoints() - 1 || j != 2) {
        contents->append(delimiter);
      }
    }
  }
}

void IGESWriter::GetParameterData2D(std::string *contents,
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
                {"128", std::to_string(spl->GetKnotVector(0)->GetNumberOfKnots() - spl->GetDegree(0).Get() - 2),
                 std::to_string(spl->GetKnotVector(1)->GetNumberOfKnots() - spl->GetDegree(1).Get() - 2),
                 std::to_string(spl->GetDegree(0).Get()), std::to_string(spl->GetDegree(1).Get()),
                 "0", "0", std::to_string(isRational), "0", "0"}, delimiter);
  auto knots1 = *spl->GetKnotVector(0);
  auto knots2 = *spl->GetKnotVector(1);
  for (auto &i : knots1) {
    contents->append(util::string_operations::GetStringWithHighPrecision(i.Get()) + delimiter);
  }
  for (auto &i : knots2) {
    contents->append(util::string_operations::GetStringWithHighPrecision(i.Get()) + delimiter);
  }
  for (int i = 0; i < spl->GetTotalNumberOfControlPoints(); ++i) {
    contents->append(util::string_operations::GetStringWithHighPrecision(spl->GetWeight(i)) + delimiter);
  }
  for (int i = 0; i < spl->GetTotalNumberOfControlPoints(); ++i) {
    for (int j = 0; j < 3; ++j) {
      contents->append(util::string_operations::GetStringWithHighPrecision(Get3DControlPoint(spl, i, j)));
      if (i != spl->GetTotalNumberOfControlPoints() - 1 || j != 2) {
        contents->append(delimiter);
      }
    }
  }
}

template<int PARAMETRIC_DIMENSIONALITY>
double IGESWriter::Get3DControlPoint(std::shared_ptr<spl::Spline<PARAMETRIC_DIMENSIONALITY>> spline, int index,
    int direction) const {
  if (direction < spline->GetPointDim()) {
    return spline->GetControlPoint(index, direction);
  }
  return 0.0;
}

std::vector<std::string> IGESWriter::GetDataEntry(int paramStart,
                                                      int paramLength,
                                                      const std::any &spline,
                                                      int *dLine) const {
  std::string contents;
  AddToContents(&contents,
                {GetBlock(std::to_string(GetDimension(spline)), 8, true), GetBlock(std::to_string(paramStart), 8, true),
                 GetBlock("0", 8, true), GetBlock("0", 8, true), GetBlock("0", 8, true), GetBlock("0", 8, true),
                 GetBlock("0", 8, true), GetBlock("0", 8, true), " 0 0 0 1",
                 GetBlock(std::to_string(GetDimension(spline)), 8, true), GetBlock("0", 8, true),
                 GetBlock("0", 8, true), GetBlock(std::to_string(paramLength), 8, true), GetBlock("0", 8, true),
                 GetBlock("", 16, true), GetBlock("SPLINE", 8, true), GetBlock("0", 8, true)}, "");
  return GetDataEntrySectionLayout(&contents, dLine);
}

std::vector<std::string> IGESWriter::GetDataEntrySectionLayout(const std::string *contents, int *dLine) const {
  std::vector<std::string> dataEntrySection;
  for (auto i = 0u; i <= (contents->size() - 1) / 72; i++) {
    dataEntrySection.emplace_back(GetBlock(contents->substr(i * 72, 72), 72, false)
                                      + 'D' + GetBlock(std::to_string(++(*dLine)), 7, true));
  }
  return dataEntrySection;
}

std::vector<std::string> IGESWriter::GetTerminateSection(int linesS, int linesG, int linesD, int linesP) const {
  return {'S' + GetBlock(std::to_string(linesS), 7, true) + 'G' + GetBlock(std::to_string(linesG), 7, true) + 'D'
              + GetBlock(std::to_string(linesD), 7, true) + 'P' + GetBlock(std::to_string(linesP), 7, true)
              + GetBlock("T", 41, true) + GetBlock("1", 7, true)};
}

int IGESWriter::GetDimension(const std::any &spline) const {
  int spline_dimension = 0;
  try {
    spline_dimension = util::any_casts::GetSplineDimension<2>(spline);
    if (spline_dimension == 1) return 126;
    return 128;
  } catch (std::logic_error &message) {
    throw std::runtime_error("Only splines of dimensions 1 or 2 can be written to an iges file.");
  }
}

bool IGESWriter::IsRational(std::any spline) const {
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

void IGESWriter::AddToContents(std::string *contents,
                                   const std::vector<std::string> &add,
                                   const std::string &delimiter) const {
  for (const auto &i : add) {
    contents->append(i + delimiter);
  }
}

std::string IGESWriter::GetHollerithFormat(const std::string &str) const {
  return std::to_string(str.size()) + "H" + str;
}

std::string IGESWriter::GetBlock(const std::string &str, int width, bool right) const {
  return right ? (std::string(width - str.size(), ' ') + str) : (str + std::string(width - str.size(), ' '));
}

void IGESWriter::WriteFile(std::ofstream &file,
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

void IGESWriter::AppendToFile(std::ofstream &file, const std::vector<std::string> &contents) const {
  for (auto &content : contents) {
    file << content << std::endl;
  }
}

std::string IGESWriter::GetTime() const {
  struct tm ptm = util::system_operations::GetTime();
  std::string date = std::to_string((ptm.tm_year + 1900) * 10000 + (ptm.tm_mon + 1) * 100 + ptm.tm_mday);
  std::string time = std::to_string(ptm.tm_hour * 10000 + ptm.tm_min * 100 + ptm.tm_sec);
  return GetHollerithFormat(date + "." + time);
}

double IGESWriter::GetHighestValue(std::vector<std::any> splines) const {
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
}  // namespace splinelib::src::io
