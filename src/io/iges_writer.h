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

#include <any>
#include <ctime>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <config.h>
#include "spline.h"
#include "b_spline.h"
#include "multi_index_handler.h"

namespace io {
class IGESWriter {
 public:
  IGESWriter() = default;

  void WriteIGESFile(spl::BSpline<1> b_spline, const std::string &filename) {
    std::ofstream newFile;
    newFile.open(filename);
    if (newFile.is_open()) {
      WriteFile(newFile, GetStartSection(), GetGlobalSection(filename, ",", ";", b_spline), GetParameterData(",", ";", b_spline));
    } else {
      throw std::runtime_error("The IGES file couldn't be opened.");
    }
    newFile.close();
  }

 private:
  std::vector<std::string> GetStartSection() {
    return {GetBlock("IGES file generated by SplineLib", 72, false) + "S" + GetBlock(GetString(1), 7, true)};
  }

  std::vector<std::string> GetGlobalSection(const std::string &filename, const std::string &delimiter,
                                            const std::string &endDelimiter, spl::BSpline<1> b_spline) {
    std::string contents;
    AddToContents(contents, {GetHollerithFormat(delimiter), GetHollerithFormat(endDelimiter),
                             GetHollerithFormat("unknown"), filename, GetHollerithFormat("SplineLib"),
                             GetHollerithFormat("pre_release"),
                             GetString(std::numeric_limits<int>::digits),
                             GetString(std::numeric_limits<float>::max_exponent),
                             GetString(std::numeric_limits<float>::digits),
                             GetString(std::numeric_limits<double>::max_exponent),
                             GetString(std::numeric_limits<double>::digits), GetHollerithFormat("unknown"),
                             GetString(1.0),
                             GetString(2), GetHollerithFormat("MM"), GetString(1), GetString(0.001),
                             GetHollerithFormat(GetTime()), GetString(0.000001), GetString(GetHighestValue(b_spline)),
                             GetHollerithFormat("unknown"), GetHollerithFormat("unknown"), GetString(11),
                             GetString(0), GetHollerithFormat(GetTime())}, delimiter);
    contents += endDelimiter;
    return GetGlobalSectionLayout(contents);
  }

  std::vector<std::string> GetGlobalSectionLayout(const std::string &contents) {
    int line = 0;
    std::vector<std::string> globalSection;
    for (unsigned long i = 0; i <= contents.size() / 72; i++) {
      globalSection.emplace_back(GetBlock(contents.substr(i * 72, 72), 72, false)
                                     + 'G' + GetBlock(GetString(++line), 7, true));
    }
    return globalSection;
  }

  std::vector<std::string> GetParameterData(const std::string &delimiter, const std::string &endDelimiter,
                                            spl::BSpline<1> b_spline) {
    std::string contents;
    GetParameterData1D(contents, delimiter, b_spline);
    contents += endDelimiter;
    return GetParameterSectionLayout(contents);
  }

  std::vector<std::string> GetParameterSectionLayout(const std::string &contents) {
    int line = 0;
    std::vector<std::string> parameterData;
    for (unsigned long i = 0; i <= contents.size() / 64; i++) {
      parameterData.emplace_back(GetBlock(contents.substr(i * 64, 64), 64, false)
                                     + ' ' + GetBlock(GetString(1), 7, true)
                                     + 'P' + GetBlock(GetString(++line), 7, true));
    }
    return parameterData;
  }


  void GetParameterData1D(std::string &contents, const std::string &delimiter, spl::BSpline<1> b_spline) {
    AddToContents(contents,
                  {GetString(126),
                   GetString(b_spline.GetKnotVector(0)->GetNumberOfKnots() - b_spline.GetDegree(0).get() - 2),
                   GetString(b_spline.GetDegree(0).get()), GetString(0), GetString(0), GetString(1), GetString(0)},
                  delimiter);

    auto knots = b_spline.GetKnots()[0];
    for (size_t i = 0; i < knots.size(); ++i) {
      contents += GetString(knots[i].get()) + delimiter;
    }
    std::vector<double> weights = b_spline.GetWeights();
    for (size_t i = 0; i < weights.size(); ++i) {
      contents += GetString(weights[i]) + delimiter;
    }
    std::vector<double> control_points = b_spline.GetControlPoints();
    for (size_t i = 0; i < control_points.size(); ++i) {
      contents += GetString(control_points[i]);
    }
  }

  void AddToContents(std::string &contents, const std::vector<std::string> &add, const std::string &delimiter) {
    for (int i = 0; i < add.size() - 1; ++i) {
      contents += add[i] + delimiter;
    }
  }

  std::string GetHollerithFormat(const std::string &str) {
    return GetString(str.size()) + "H" + str;
  }

  std::string GetBlock(const std::string &str, int width, bool right) {
    return right ? (std::string(width - str.size(), ' ') + str) : (str + std::string(width - str.size(), ' '));
  }

  template<typename T>
  std::string GetString(T value) {
    return std::to_string(value);
  }

  void WriteFile(std::ofstream &file, const std::vector<std::string> &start,
                 const std::vector<std::string> &global/*, const std::vector<std::string> &data*/,
                 const std::vector<std::string> &parameter/*, const std::vector<std::string> &terminate*/) {
    AppendToFile(file, start);
    AppendToFile(file, global);
    //AppendToFile(file, data);
    AppendToFile(file, parameter);
    //AppendToFile(file, terminate);
  }

  void AppendToFile(std::ofstream &file, const std::vector<std::string> &contents) {
    for (auto &content : contents) {
      file << content << std::endl;
    }
  }

  std::string GetTime() {
    time_t timer;
    time(&timer);
    struct tm *ptm;
    ptm = localtime(&timer);
    std::string
        date = GetString((ptm->tm_year + 1900) * 10000 + (ptm->tm_mon + 1) * 100 + ptm->tm_mday);
    std::string time = GetString(ptm->tm_hour * 10000 + ptm->tm_min * 100 + ptm->tm_sec);
    if (time.size() == 5) {
      time = '0' + time;
    }
    return date;
  }

  double GetHighestValue(spl::BSpline<1> b_spline) {
    return b_spline.GetExpansion();
  }
};
}

#endif  // SRC_IO_IGES_WRITER_H_
