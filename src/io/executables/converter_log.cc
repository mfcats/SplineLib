/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "converter_log.h"

#include <ctime>
#include <fstream>
#include <iostream>

#include "string_operations.h"
#include "system_operations.h"
#include "vector_utils.h"

namespace splinelib::src::io {
ConverterLog::ConverterLog() : log_file_("") {}

ConverterLog::ConverterLog(const char *log_file) : log_file_(log_file) {
  std::ifstream log;
  log.open(log_file);
  if (!log.good()) {
    throw std::runtime_error("Log file could not be opened.");
  }
  std::string line;
  while (getline(log, line)) {
    if (util::StringOperations::StartsWith(line, "input:")) {
      getline(log, input_);
    }
    if (util::StringOperations::StartsWith(line, "output:")) {
      getline(log, output_);
    }
    if (util::StringOperations::StartsWith(line, "options:")) {
      while (getline(log, line) && !util::StringOperations::EndsWith(util::StringOperations::trim(line), ":")) {
        if (line == "all") {
          positions_.clear();
          positions_.emplace_back(-1);
          break;
        }
        positions_.push_back(std::stoi(line));
      }
    }
    if (util::StringOperations::StartsWith(line, "scattering:")) {
      while (getline(log, line) && !util::StringOperations::EndsWith(util::StringOperations::trim(line), ":")) {
        std::vector<std::string> strings = util::StringOperations::split(line, ' ');
        scattering_.emplace_back(util::StringOperations::StringVectorToNumberVector<int>(strings));
      }
    }
  }
}

const char *ConverterLog::GetInput() const {
  return input_.c_str();
}

const char *ConverterLog::GetOutput() const {
  return output_.c_str();
}

std::vector<int> ConverterLog::GetPositions(std::vector<int> possible_positions) {
  if (positions_[0] == -1) {
    written_ = possible_positions;
  } else {
    for (const auto &pos : positions_) {
      if (std::find(possible_positions.begin(), possible_positions.end(), pos) != possible_positions.end()) {
        written_.emplace_back(pos);
      } else {
        not_written_.emplace_back(pos);
      }
    }
  }
  return written_;
}

std::vector<std::vector<int>> ConverterLog::GetScattering() {
  std::vector<int> indices;
  if (positions_[0] != -1) {
    for (int i = 0; i < static_cast<int>(positions_.size()); ++i) {
      if (std::find(written_.begin(), written_.end(), positions_[i]) != written_.end()) {
        indices.emplace_back(i);
      }
    }
  } else {
    indices = written_;
  }
  return util::VectorUtils<std::vector<int>>::FilterVector(scattering_, indices);
}

void ConverterLog::WriteLog() const {
  std::ofstream log;
  log.open(log_file_, std::ios_base::app);
  log << "\n\nlog:\n" << GetTime();
  if (!written_.empty()) {
    log << "The " << (OneSpline() ? "spline" : "splines") << " at position" << (OneSpline() ? " " : "s ");
    for (auto i = 0u; i < written_.size() - 1; ++i) {
      log << written_[i] << ", ";
    }
    log << written_.back();
  } else {
    log << "No splines";
  }
  log << " in file " << input_ << (OneSpline() ? " has" : " have") << " been written to " << output_;
}

void ConverterLog::PrintHelp() const {
  std::cout << "The log file has to be of the following format:" << std::endl
            << "input:\n# path to input file\n" << std::endl
            << "output:\n# path to output file\n" << std::endl
            << "options:\n# list separated by line breaks of positions of splines in input file to be written to the "
            << "output file, 'all' for all splines in input file\n" << std::endl
            << "If this is a converter to VTK format, the log file has to expanded by the following entry:" << std::endl
            << "scattering:\n# scattering for each spline separated by line breaks and for each dimension separated "
            << "by spaces\n" << std::endl
            << "The log file is expanded by a log entry if converting succeed:" << std::endl
            << "log:\n# time date\n# spline positions in input file that have been written to output file.\n";
}

std::string ConverterLog::GetTime() const {
  struct tm timeinfo = util::SystemOperations::GetTime();
  return asctime(&timeinfo);
}

bool ConverterLog::OneSpline() const {
  return written_.size() == 1;
}
}  // namespace splinelib::src::io
