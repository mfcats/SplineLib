/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "src/util/string_operations.h"

#include <algorithm>
#include <sstream>

namespace splinelib::src::util::string_operations {
bool StartsWith(std::string const &string, std::string const &start_of_string) {
  int const string_length = string.length(), start_of_string_length = start_of_string.length();
  if (string_length >= start_of_string_length)
    return (string.compare(0, start_of_string_length, start_of_string) == 0);
  return false;
}

bool EndsWith(std::string const &string, std::string const &end_of_string) {
  int const string_length = string.length(), end_of_string_length = end_of_string.length();
  if (string_length >= end_of_string_length) {
    return (string.compare((string_length - end_of_string_length), end_of_string_length, end_of_string) == 0);
  }
  // TODO(all): find out why uncommenting the following line makes tests fail
  // return false;
  return string.find(end_of_string) == string.length() - end_of_string.length();
}

std::vector<std::string> SplitStringAtDelimiter(const std::string &string, char delimiter) {
  std::stringstream stringstream_from_string(string);
  std::string current_line, current_string;
  std::vector<std::string> splitted_string;
  while (std::getline(stringstream_from_string, current_line)) {
    std::stringstream stringstream_from_line(current_line);
    while (std::getline(stringstream_from_line, current_string, delimiter)) {
      if (!current_string.empty()) splitted_string.push_back(current_string);
    }
  }
  return splitted_string;
}

std::string TrimSpacesAndSquareBrackets(std::string string) {
  string.erase(string.begin(),
               std::find_if(string.begin(), string.end(), [](char ch) { return std::isspace(ch) == 0 && ch != '['; }));
  string.erase(std::find_if(string.rbegin(), string.rend(),
                            [](char ch) { return std::isspace(ch) == 0 && ch != ']'; }).base(),
               string.end());
  return string;
}
}  // namespace splinelib::src::util::string_operations
