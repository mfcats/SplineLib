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
  static std::vector<std::string> split(const std::string &string, char delimiter) {
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

  static bool StartsWith(const std::string &string, const std::string &start_of_string) {
    return string.find(start_of_string) == 0;
  }

  static bool EndsWith(const std::string &string, const std::string &end_of_string) {
    return string.find(end_of_string) == string.length() - end_of_string.length();
  }

  static double StringToDouble(std::string string) {
    std::vector<int> left;
    std::vector<int> right;
    int sign = 1;
    if (string[0] == '-') {
      sign = -1;
      string.erase(0, 1);
    }
    std::size_t found = string.find_first_of('.');
    for (auto i = 0u; i < string.length(); ++i) {
      if (i < found) {
        left.emplace_back(std::stoi(string.substr(i, 1)));
      } else if (i > found) {
        right.emplace_back(std::stoi(string.substr(i, 1)));
      }
    }
    double result = 0;
    double factor = 1 * (std::pow(10, left.size() - 1));
    for (int i : left) {
      result += i * factor;
      factor /= 10;
    }
    factor = 0.1;
    for (int i : right) {
      result += i * factor;
      factor /= 10;
    }
    return sign * result;
  }
};
}  // namespace util

#endif  // SRC_UTIL_STRING_OPERATIONS_H_
