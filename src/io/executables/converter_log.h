/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/
#ifndef SRC_IO_EXECUTABLES_CONVERTER_LOG_H_
#define SRC_IO_EXECUTABLES_CONVERTER_LOG_H_

#include <string>
#include <vector>

namespace splinelib::src::io {
class ConverterLog {
 public:
  ConverterLog();
  explicit ConverterLog(const char *log_file);

  const char *GetInput() const;
  const char *GetOutput() const;
  std::vector<int> GetPositions(std::vector<int> possible_positions);
  std::vector<std::vector<int>> GetScattering();

  void WriteLog() const;
  void PrintHelp() const;

 private:
  std::string GetTime() const;
  bool OneSpline() const;

  const char *log_file_;
  std::string input_;
  std::string output_;
  std::vector<int> positions_;
  std::vector<int> written_;
  std::vector<int> not_written_;
  std::vector<std::vector<int>> scattering_;
};
}  // namespace splinelib::src::splinelib::src::io

#endif  // SRC_IO_EXECUTABLES_CONVERTER_LOG_H_
