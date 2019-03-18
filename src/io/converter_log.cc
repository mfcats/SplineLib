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

#include "string_operations.h"
#include "vector_utils.h"

io::ConverterLog::ConverterLog(const char *log_file) : log_file_(log_file) {
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

const char *io::ConverterLog::GetInput() const {
  return input_.c_str();
}

const char *io::ConverterLog::GetOutput() const {
  return output_.c_str();
}

std::vector<int> io::ConverterLog::GetPositions(std::vector<int> possible_positions) {
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

std::vector<std::vector<int>> io::ConverterLog::GetScattering() {
  return util::VectorUtils<std::vector<int>>::FilterVector(scattering_, written_);
}

void io::ConverterLog::WriteLog() {
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

std::string io::ConverterLog::GetTime() {
  struct tm timeinfo{};
  time_t rawtime;
  rawtime = time(&rawtime);
  localtime_r(&rawtime, &timeinfo);
  return asctime(&timeinfo);
}

bool io::ConverterLog::OneSpline() {
  return written_.size() == 1;
}
