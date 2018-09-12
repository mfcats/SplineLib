/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_UTIL_STRING_OPERATIONS_H_
#define SRC_UTIL_STRING_OPERATIONS_H_

#include <string>
#include <vector>

namespace util {
class StringOperations {
 public:
  static
  std::vector<std::string> split(const std::string &string, char delimiter) {
    std::stringstream ss_string(string);
    std::string current_line;
    std::string current_value;
    std::vector<std::string> splitted_string;
    while (std::getline(ss_string, current_line)) {
      std::stringstream ss_current_line(current_line);
      while (std::getline(ss_current_line, current_value, delimiter)) {
        if (!current_value.empty()) splitted_string.push_back(current_value);
      }
    }
    return splitted_string;
  }

  template<class T>
  static std::vector<T> StringVectorToNumberVector(const std::vector<std::string> &string_vector) {
    std::vector<T> converted;
    for (const std::string &string : string_vector) {
      converted.emplace_back(strtod(string.c_str(), nullptr));
    }
    return converted;
  }
};
}  // namespace util

#endif  // SRC_UTIL_STRING_OPERATIONS_H_
