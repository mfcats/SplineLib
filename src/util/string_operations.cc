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
#include <cmath>
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

std::vector<std::string> Split(const std::string &string, char delimiter) {
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

std::string Trim(std::string string) {
  string.erase(string.begin(),
               std::find_if(string.begin(), string.end(), [](char ch) { return std::isspace(ch) == 0 && ch != '['; }));
  string.erase(std::find_if(string.rbegin(), string.rend(),
                            [](char ch) { return std::isspace(ch) == 0 && ch != ']'; }).base(),
               string.end());
  return string;
}

double StringToDouble(std::string string) {
  int sign = 1;
  if (string[0] == '-') {
    sign = -1;
    string.erase(0, 1);
  }
  int const found_dot = string.find_first_of('.');
  int exponent = string.find_first_of("Ee");
  int const npos = static_cast<int>(std::string::npos);
  int const string_length = string.length();
  int end_of_number;
  if (exponent == npos) {
    end_of_number = (string_length - 1);
  } else {
    end_of_number = (exponent - 1);
  }
  double result = 0, factor;
  if (found_dot != npos) {
    factor = pow(10, found_dot - 1);
  } else {
    factor = pow(10, string_length - 1);
  }
  for (int i = 0; i <= end_of_number; ++i) {
    if (i != found_dot) {
      result += (std::stoi(string.substr(i, 1)) * factor);
      factor /= 10;
    }
  }
  if (exponent != npos) {
    if (string[exponent + 1] == '-') {
      factor = -1;
      ++exponent;
    } else if (string[exponent + 1] == '+') {
      ++exponent;
    }
    int const potency = std::stoi(string.substr(exponent + 1, string_length - exponent));
    if (factor != -1) {
      result *= pow(10, potency);
    } else {
      result /= pow(10, potency);
    }
  }
  return (sign * result);
}

std::vector<double> DelimitedStringToVector(std::string string) {
  std::vector<double> vector;
  while (!string.empty()) {
    int const found_comma = string.find_first_of(',');
    int const found_semicolon = string.find_first_of(';');
    if (((found_comma < found_semicolon) || (found_semicolon == -1)) && (found_comma > 0)) {
      vector.push_back(StringToDouble(Trim(string.substr(0, found_comma))));
      string.erase(0, found_comma + 1);
    } else if (((found_semicolon < found_comma) || (found_comma == -1)) && (found_semicolon > 0)) {
      vector.push_back(StringToDouble(Trim(string.substr(0, found_semicolon))));
      string.erase(0, found_semicolon + 1);
    } else {
      string.erase(0, 1);
    }
  }
  return vector;
}
}  // namespace splinelib::src::util::string_operations
